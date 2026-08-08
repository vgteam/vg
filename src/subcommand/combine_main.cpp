/** \file combine_main.cpp
 *
 * Defines the "vg combine" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include <vg/io/vpkg.hpp>
#include <vg/io/message_iterator.hpp>
#include <vg/io/blocked_gzip_input_stream.hpp>

#include "subcommand.hpp"

#include "../handle.hpp"
#include "../vg.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_combine(char** argv) {
    cerr << "usage: " << argv[0] << " combine [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
         << "Combines one or more graphs into a single file, regardless of input format." << endl
         << "Node IDs will be modified as needed to resolve conflicts (same as vg ids -j)." << endl
         << endl
         << "Options:" << endl
         << "  -c, --cat-proto       merge by converting each to Protobuf (if not already)" << endl
         << "                        and catting the results." << endl
         << "                        node IDs not modified [DEPRECATED]" << endl
         << "  -p, --connect-paths   add edges necessary to connect paths with the same name" << endl
         << "                        which are present in different graphs." << endl
         << "                        ex: If path x is present in graphs N-1 and N, then" << endl
         << "                        an edge connecting the last node of x in N-1 " << endl
         << "                        and the first node of x in N will be added." << endl
         << "  -f, --connect-fragments" << endl
         << "                        [DEFAULT MODE] treat paths sharing" << endl
         << "                        sample/locus/haplotype/phase block but differing only by" << endl
         << "                        subrange (e.g. CHM13#0#chr22 and" << endl
         << "                        CHM13#0#chr22[17475777]) as fragments of one path. Sort" << endl
         << "                        fragments by start offset, concatenate in order, add" << endl
         << "                        the missing inter-fragment edges, and rewrite the path" << endl
         << "                        name to a single [start-end] subrange covering the" << endl
         << "                        union. This runs by default whenever none of -c/-p/-m" << endl
         << "                        is given." << endl
         << "  -s, --shared-nodes    modify the default (-f) mode: do NOT renumber" << endl
         << "                        IDs when merging. Node IDs that appear in both" << endl
         << "                        inputs must carry the same sequence; one copy is" << endl
         << "                        kept (intended for `vg chunk` output that shares" << endl
         << "                        boundary nodes across chunks). Fragments are merged" << endl
         << "                        into one path, trimming duplicated boundary steps." << endl
         << "                        Valid only with the default mode; cannot be" << endl
         << "                        combined with -c/-p/-m." << endl
         << "  -m, --merge           combine chunked-assembly graphs that share" << endl
         << "                        REFERENCE-sense paths (e.g. CHM13#0#chr22)" << endl
         << "                        and overlap on them. Buckets inputs by" << endl
         << "                        REFERENCE identity (sample/locus/haplotype)," << endl
         << "                        so chunks of different chromosomes (e.g." << endl
         << "                        CHM13#0#chr21 and CHM13#0#chr22) can be passed" << endl
         << "                        together; each reference is stitched on its" << endl
         << "                        own and kept as a separate path in the output." << endl
         << "                        Within a reference, sorts inputs by REFERENCE" << endl
         << "                        start offset, validates that each chunk's" << endl
         << "                        boundary is a single chain node shared by" << endl
         << "                        every path (and equal to the REFERENCE end)," << endl
         << "                        trims each right chunk's leading node by the" << endl
         << "                        REFERENCE-offset overlap, connects the left" << endl
         << "                        end to the trimmed right start, then merges" << endl
         << "                        fragments like -f." << endl
         << "  -u, --fusion          requires -m. Instead of just adding an edge" << endl
         << "                        between the left end node and the right's" << endl
         << "                        trimmed start node, fuse them into a single" << endl
         << "                        node carrying the concatenated sequence." << endl
         << "  -h, --help            print this help message to stderr and exit" << endl;
}

static int cat_proto_graphs(int argc, char** argv, const Logger& logger);
static int connect_paths_combine(int argc, char** argv, const Logger& logger);
static int connect_fragments_combine(int argc, char** argv, const Logger& logger,
                                     bool shared_nodes);
static int merge_combine(int argc, char** argv, const Logger& logger, bool fusion);

int main_combine(int argc, char** argv) {
    Logger logger("vg combine");

    if (argc == 2) {
        help_combine(argv);
        return 1;
    }

    bool connect_paths = false;
    bool connect_fragments = false;
    bool shared_nodes = false;
    bool cat_proto = false;
    bool merge = false;
    bool fusion = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"connect-paths", no_argument, 0, 'p'},
            {"connect-fragments", no_argument, 0, 'f'},
            {"shared-nodes", no_argument, 0, 's'},
            {"cat-proto", no_argument, 0, 'c'},
            {"merge", no_argument, 0, 'm'},
            {"fusion", no_argument, 0, 'u'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?pfscmu",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            connect_paths = true;
            break;
        case 'f':
            connect_fragments = true;
            break;
        case 's':
            shared_nodes = true;
            break;
        case 'c':
            cat_proto = true;
            break;
        case 'm':
            merge = true;
            break;
        case 'u':
            fusion = true;
            break;
        case 'h':
        case '?':
            help_combine(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (fusion && !merge) {
        logger.error() << "--fusion/-u requires --merge/-m." << endl;
    }

    if (cat_proto) {
        if (connect_paths)
        logger.warn() << "--cat-proto/-c option is deprecated "
                      << "and will be removed in a future version of vg." << endl;
        return cat_proto_graphs(argc, argv, logger);
    }

    if (merge) {
        if (connect_paths || connect_fragments || shared_nodes) {
            logger.error() << "--merge/-m is mutually exclusive with --connect-paths/-p, "
                           << "--connect-fragments/-f, and --shared-nodes/-s." << endl;
        }
        return merge_combine(argc, argv, logger, fusion);
    }

    if (connect_paths) {
        // --connect-paths/-p mode (explicit), mutually exclusive with -f and -s.
        if (connect_fragments || shared_nodes) {
            logger.error() << "--connect-paths/-p is mutually exclusive with "
                           << "--connect-fragments/-f and --shared-nodes/-s." << endl;
        }
        return connect_paths_combine(argc, argv, logger);
    }

    // --connect-fragments/-f is the DEFAULT mode. Unless the user explicitly
    // selects another mode (--cat-proto/-c and --merge/-m are handled above,
    // --connect-paths/-p is handled just above), fragment merging runs. Passing
    // -f explicitly is equivalent to the default.
    //
    // --shared-nodes/-s is a modifier of this default mode: it is forwarded so
    // the per-graph copy keeps source IDs intact and shared boundary nodes are
    // deduplicated before fragments are merged into a single coordinate-spanning
    // path. It is therefore only valid here (not with -c/-p/-m).
    return connect_fragments_combine(argc, argv, logger, shared_nodes);
}

// Register subcommand
static Subcommand vg_combine("combine", "merge multiple graph files together", main_combine);

// --connect-paths / -p mode.
//
// Loads the first graph as the destination accumulator, then for each
// subsequent graph renumbers its node IDs out of the way and appends it,
// connecting same-named paths end-to-end across graphs (an edge is added from
// the last node of a path in the accumulator to the first node of the same-named
// path in the incoming graph). See handlealgs::append_path_handle_graph.
int connect_paths_combine(int argc, char** argv, const Logger& logger) {

    unique_ptr<MutablePathMutableHandleGraph> first_graph;
    string first_graph_filename = get_input_file_name(optind, argc, argv);
    first_graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(first_graph_filename);
    int64_t max_node_id = first_graph->max_node_id();

    while (optind < argc) {

        unique_ptr<MutablePathMutableHandleGraph> graph;
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_filename);

        // join the id spaces if necessary
        int64_t delta = max_node_id - graph->min_node_id();
        if (delta >= 0) {
            graph->increment_node_ids(delta + 1);
        }
        max_node_id = graph->max_node_id();

        // -p: append, connecting same-named paths end-to-end across graphs.
        handlealgs::append_path_handle_graph(graph.get(), first_graph.get(), true);
    }

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(first_graph.get(), cout);

    return 0;
}

// Helper: compute the total sequence length spanned by a path in `graph`.
static size_t compute_path_length(const PathHandleGraph* graph, const path_handle_t& path) {
    size_t length = 0;
    for (handle_t h : graph->scan_path(path)) {
        length += graph->get_length(h);
    }
    return length;
}

// Bucket paths in `dest` by their non-coordinate metadata
// (sense, sample, locus, haplotype, circular) -- phase_block is intentionally
// excluded so paths that use phase_block as a start coordinate land in the
// same bucket -- then sort each bucket by start offset, drop the original
// fragment paths, and re-create a single coordinate-spanning path per group.
// Inter-fragment edges are added where missing; runs of duplicated handles
// at the seam (suffix-of-prev == prefix-of-next) are trimmed so the merged
// path doesn't revisit shared boundary nodes.
//
// This is the second pass of -f / -s+-f / -m: it expects `dest` to already
// contain every input graph's nodes, edges, and per-fragment paths.
static void merge_path_fragments(MutablePathMutableHandleGraph* dest,
                                 const Logger& logger) {
    using handlegraph::PathSense;
    using handlegraph::PathMetadata;
    using handlegraph::subrange_t;
    using handlegraph::offset_t;

    struct Fragment {
        path_handle_t handle;
        string name;
        subrange_t subrange;
        size_t phase_block;
        size_t length;
    };
    using GroupKey = std::tuple<PathSense, string, string, size_t, bool>;
    std::map<GroupKey, std::vector<Fragment>> groups;

    dest->for_each_path_handle([&](const path_handle_t& path_handle) {
        GroupKey key{
            dest->get_sense(path_handle),
            dest->get_sample_name(path_handle),
            dest->get_locus_name(path_handle),
            dest->get_haplotype(path_handle),
            dest->get_is_circular(path_handle)
        };
        Fragment f;
        f.handle = path_handle;
        f.name = dest->get_path_name(path_handle);
        f.subrange = dest->get_subrange(path_handle);
        f.phase_block = dest->get_phase_block(path_handle);
        f.length = compute_path_length(dest, path_handle);
        groups[key].push_back(std::move(f));
    });

    auto fragment_start = [](const Fragment& f) -> offset_t {
        if (f.subrange != PathMetadata::NO_SUBRANGE) {
            return f.subrange.first;
        }
        if (f.phase_block != PathMetadata::NO_PHASE_BLOCK) {
            return (offset_t)f.phase_block;
        }
        return 0;
    };
    auto fragment_end = [&](const Fragment& f) -> offset_t {
        if (f.subrange != PathMetadata::NO_SUBRANGE &&
            f.subrange.second != PathMetadata::NO_END_POSITION) {
            return f.subrange.second;
        }
        return fragment_start(f) + (offset_t)f.length;
    };

    for (auto& kv : groups) {
        auto& frags = kv.second;
        if (frags.size() <= 1) {
            continue;
        }
        const auto& key = kv.first;
        PathSense input_sense  = std::get<0>(key);
        const string& sample   = std::get<1>(key);
        const string& locus    = std::get<2>(key);
        size_t haplotype       = std::get<3>(key);
        bool is_circular       = std::get<4>(key);

        std::sort(frags.begin(), frags.end(),
                  [&](const Fragment& a, const Fragment& b) {
                      return fragment_start(a) < fragment_start(b);
                  });

        offset_t combined_start = fragment_start(frags.front());
        offset_t combined_end = 0;
        for (const auto& f : frags) {
            combined_end = std::max(combined_end, fragment_end(f));
        }
        subrange_t combined_subrange{combined_start, combined_end};

        // Method A: preserve the input sense. A HAPLOTYPE path cannot exist
        // without a phase_block (create_path_name enforces this), so instead of
        // demoting it to REFERENCE/GENERIC we keep it HAPLOTYPE and attach a
        // placeholder phase_block of 0. The real coordinates live in the
        // subrange, yielding names like HG002#2#chr22#0[26540133-27549291].
        // REFERENCE/GENERIC paths must NOT carry a phase_block, so they keep
        // NO_PHASE_BLOCK as before.
        PathSense merged_sense = input_sense;
        size_t merged_phase_block = (input_sense == PathSense::HAPLOTYPE)
            ? 0
            : PathMetadata::NO_PHASE_BLOCK;

        string merged_name = PathMetadata::create_path_name(
            merged_sense, sample, locus, haplotype,
            merged_phase_block, combined_subrange);

        if (dest->has_path(merged_name)) {
            bool is_own_fragment = false;
            for (const auto& f : frags) {
                if (f.name == merged_name) { is_own_fragment = true; break; }
            }
            if (!is_own_fragment) {
                logger.error() << "Merged fragment name \"" << merged_name
                               << "\" collides with an existing unrelated path in the combined graph."
                               << endl;
            }
        }

        std::vector<std::vector<handle_t>> fragment_steps;
        fragment_steps.reserve(frags.size());
        for (const auto& f : frags) {
            std::vector<handle_t> steps;
            for (handle_t h : dest->scan_path(f.handle)) {
                steps.push_back(h);
            }
            fragment_steps.push_back(std::move(steps));
        }

        for (const auto& f : frags) {
            dest->destroy_path(f.handle);
        }

        path_handle_t merged = dest->create_path(
            merged_sense, sample, locus, haplotype,
            merged_phase_block, combined_subrange, is_circular);

        bool have_prev = false;
        handle_t prev_tail{};
        const std::vector<handle_t>* prev_steps_ptr = nullptr;
        for (size_t i = 0; i < fragment_steps.size(); ++i) {
            const auto& steps = fragment_steps[i];
            if (steps.empty()) {
                continue;
            }
            size_t skip = 0;
            if (have_prev && prev_steps_ptr != nullptr) {
                const auto& prev_steps = *prev_steps_ptr;
                size_t max_overlap = std::min(prev_steps.size(), steps.size());
                for (size_t k = max_overlap; k >= 1; --k) {
                    bool match = true;
                    for (size_t off = 0; off < k; ++off) {
                        if (prev_steps[prev_steps.size() - k + off] != steps[off]) {
                            match = false;
                            break;
                        }
                    }
                    if (match) {
                        skip = k;
                        break;
                    }
                }
            }
            if (have_prev && skip == 0 && !dest->has_edge(prev_tail, steps.front())) {
                dest->create_edge(prev_tail, steps.front());
            }
            for (size_t j = skip; j < steps.size(); ++j) {
                dest->append_step(merged, steps[j]);
            }
            prev_tail = steps.back();
            prev_steps_ptr = &steps;
            have_prev = true;
        }

        logger.info() << "Merged " << frags.size() << " fragments into \""
                      << merged_name << "\"" << endl;
    }
}

// --connect-fragments / -f mode.
//
// Combines multiple graphs and merges paths that are fragments of the same
// logical path. Paths are bucketed by (sense, sample, locus, haplotype,
// circular). phase_block is *not* part of the key: when it varies inside a
// bucket it is interpreted as a start coordinate (the PathMetadata spec
// explicitly allows phase_block to be "a number or a start coordinate"), so
// e.g. `HG002#1#chr22#0` and `HG002#1#chr22#15419374` are recognized as
// fragments of the same logical chr22 and merged.
//
// Each fragment's effective coordinate range on the merged path is derived as:
//   - No subrange, no phase_block -> [0, path-length)
//   - No subrange, phase_block=s  -> [s, s + path-length)        (single offset)
//   - subrange=[s]                -> [s, s + path-length)        (single offset)
//   - subrange=[s-e]              -> [s, e)                      (range)
// Fragments are then sorted by start, the missing inter-fragment edges are
// added, and the union [min start, max end-exclusive) becomes the merged path's
// subrange. The merged path drops phase_block entirely; if the inputs were
// HAPLOTYPE-sense purely because of phase_block, the merged path becomes
// REFERENCE-sense (or GENERIC if there is no sample).
//
// When `shared_nodes` is true (-s + -f), the per-graph copy preserves source
// node IDs and deduplicates nodes that already exist in the accumulator
// (sequence-verified). The fragment-stitching loop also checks for shared
// boundary handles between adjacent fragments and trims the duplicate steps,
// so the merged path traverses each shared node only once.
int connect_fragments_combine(int argc, char** argv, const Logger& logger,
                              bool shared_nodes) {

    using handlegraph::nid_t;

    // (1) Load the first graph; it doubles as the destination accumulator.
    string first_graph_filename = get_input_file_name(optind, argc, argv);
    unique_ptr<MutablePathMutableHandleGraph> dest =
        vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(first_graph_filename);
    int64_t max_node_id = dest->max_node_id();

    // (2) For each subsequent graph: copy nodes/edges into the destination,
    //     then copy each path with its original (per-fragment) name. The
    //     fragment merging happens in a second pass.
    //
    //     Default behavior renumbers the incoming graph's IDs out of the way.
    //     Under -s (shared_nodes) we keep source IDs and dedupe colliding
    //     nodes by ID after verifying their sequences match.
    while (optind < argc) {
        unique_ptr<MutablePathMutableHandleGraph> graph;
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_filename);

        if (shared_nodes) {
            graph->for_each_handle([&](const handle_t& src_h) {
                nid_t id = graph->get_id(src_h);
                string src_seq = graph->get_sequence(src_h);
                if (dest->has_node(id)) {
                    if (dest->get_sequence(dest->get_handle(id)) != src_seq) {
                        logger.error() << "Node " << id << " has conflicting sequences across "
                                       << "inputs; cannot combine in --shared-nodes mode." << endl;
                    }
                    // else: shared node already in dest; skip creation.
                } else {
                    dest->create_handle(src_seq, id);
                }
            });
            graph->for_each_edge([&](const edge_t& src_edge) {
                handle_t a = dest->get_handle(graph->get_id(src_edge.first),
                                              graph->get_is_reverse(src_edge.first));
                handle_t b = dest->get_handle(graph->get_id(src_edge.second),
                                              graph->get_is_reverse(src_edge.second));
                if (!dest->has_edge(a, b)) {
                    dest->create_edge(a, b);
                }
            });
            max_node_id = std::max(max_node_id, (int64_t)graph->max_node_id());
        } else {
            int64_t delta = max_node_id - graph->min_node_id();
            if (delta >= 0) {
                graph->increment_node_ids(delta + 1);
            }
            max_node_id = graph->max_node_id();

            // Copy nodes + edges only -- paths handled below so we can detect
            // name collisions explicitly.
            handlealgs::copy_handle_graph(graph.get(), dest.get());
        }

        graph->for_each_path_handle([&](const path_handle_t& src_path) {
            string path_name = graph->get_path_name(src_path);
            if (dest->has_path(path_name)) {
                logger.error() << "Paths with name \"" << path_name
                               << "\" found in multiple input graphs. Two fragments cannot share "
                               << "the same subrange suffix; please disambiguate the inputs." << endl;
            }
            if (shared_nodes) {
                // IDs match between graph and dest; replay steps directly.
                path_handle_t new_path = dest->create_path(
                    graph->get_sense(src_path),
                    graph->get_sample_name(src_path),
                    graph->get_locus_name(src_path),
                    graph->get_haplotype(src_path),
                    graph->get_phase_block(src_path),
                    graph->get_subrange(src_path),
                    graph->get_is_circular(src_path)
                );
                for (handle_t src_step : graph->scan_path(src_path)) {
                    handle_t dest_step = dest->get_handle(graph->get_id(src_step),
                                                         graph->get_is_reverse(src_step));
                    dest->append_step(new_path, dest_step);
                }
            } else {
                handlealgs::copy_path(graph.get(), src_path, dest.get());
            }
        });
    }

    // (3) Bucket and merge fragments (shared with -m).
    merge_path_fragments(dest.get(), logger);

    // (4) Serialize the merged graph.
    vg::io::save_handle_graph(dest.get(), cout);
    return 0;
}


// This is the original vg combine logic, which itself mimics using "cat" to join up protobuf files
// Since it relies on the Protobuf format itself, particular the ability to stream together chunks that
// would otherwise be invalid individually, it is probably never going to be ported to the handle graph
// api, which is why it's been relegated to the deprecated bin
int cat_proto_graphs(int argc, char** argv, const Logger& logger) {
    
    while (optind < argc) {
        get_input_file(optind, argc, argv, [&](istream& in) {
            // We're producing output in uncompressed, "VG"-type-tagged, VPKG Protobuf format.
            // We will check if this file is uncompressed or compressed VG-type-tagged data.
            
            if (vg::io::BlockedGzipInputStream::SmellsLikeGzip(in)) {
                // It is compressed.
                
                // Save our start position
                auto start = in.tellg();
                
                {
                    // Try decompressing.
                    vg::io::BlockedGzipInputStream decompressed(in);
                    if (decompressed.IsBGZF() && vg::io::MessageIterator::sniff_tag(decompressed) == "VG") {
                        // We have Blocked GZIP which we can potentially just forward.
                        // It looks like compressed VG Protobuf data.
                        
                        // Decompress it all to stdout, using the ZeroCopyInputStream API.
                        char* buffer = nullptr;
                        int bytes = 0;
                        while (cout && decompressed.Next((const void**) &buffer, &bytes)) {
                            // Each time we get bytes, write them to stdout.
                            cout.write(buffer, bytes);
                        }
                        
                        if (!cout) {
                            logger.error() << "Could not write decompressed data to output stream." << endl;
                        }
                        
                        // Do the next input file
                        return;
                    }
                }
                
                // We may have hit EOF.
                in.clear();
                
                // If we get here, it wasn't compressed VG Protobuf.
                // So we need to go back to the start of the file, since the decompressor read some.
                in.seekg(start);
                
            } else if (vg::io::MessageIterator::sniff_tag(in) == "VG") {
                // It isn't compressed, but it looks like uncompressed VG Protobuf.
                // Send the uncompressed data to stdout.
                cout << in.rdbuf();
                
                if (!cout) {
                    logger.error() << "Could not write raw data to output stream." << endl;
                }
                
                // Do the next input file
                return;
            }
            
            // If we get here, it isn't compressed or uncompressed VG protobuf.
            // Read it as a PathHandleGraph
            unique_ptr<PathHandleGraph> graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
            
            // Convert to vg::VG
            VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
            if (vg_graph == nullptr) {
                vg_graph = new vg::VG();
                handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
                // Give the unique_ptr ownership and delete the graph we loaded.
                graph.reset(vg_graph);
                // Make sure the paths are all synced up
                vg_graph->paths.to_graph(vg_graph->graph);
            }
            
            {
                // Save to stdout, uncompressed
                vg::io::ProtobufEmitter<Graph> emitter(cout, false);
                vg_graph->serialize_to_emitter(emitter);
                // Make sure the emitter goes away and writes before we check on the stream.
            }
            
            if (!cout) {
                logger.error() << "Could not write converted graph to output stream." << endl;
            }

        });
    }
    return 0;
}


// Per-input REFERENCE metadata for --merge/-m, shared between the bucket
// stitcher and the top-level driver.
struct MergeInput {
    string filename;
    unique_ptr<MutablePathDeletableHandleGraph> graph;
    path_handle_t ref_path;
    string ref_name;
    handle_t ref_start_handle;
    handle_t ref_end_handle;
    handlegraph::offset_t ref_start_offset = 0;
    handlegraph::offset_t ref_end_offset = 0;
    string sample;
    string locus;
    size_t haplotype;
    bool is_circular = false;
};

// Stitch one bucket of inputs that all share the same REFERENCE identity
// (sample, locus, haplotype, circular) -- i.e. chunks of a single chromosome --
// into one accumulator graph. This is steps (2)-(4) of --merge/-m:
//
//   (2) Sort inputs by REFERENCE start offset.
//   (3) Validate the chain boundary on every input:
//         - every non-first input's paths all start at its REF start node
//         - every non-last input's paths all end at its REF end node
//       Boundary nodes must be visited in forward orientation.
//   (4) Walk pairs left->right. For each right chunk:
//         - overlap = dest_REF_end_offset - src_REF_start_offset
//           (error on negative gap or if overlap >= start node length)
//         - if overlap > 0: divide_handle(start, overlap), pop_front_step on
//           every path in the right chunk, destroy the head piece. The new
//           start node is the tail piece.
//         - renumber the right chunk's IDs out of dest's range, copy
//           nodes/edges/paths into dest under their per-fragment names.
//         - connect: if --fusion, glue the left end node and the new right
//           start node into a single node with concatenated sequence;
//           otherwise just create_edge(left_end, new_right_start).
//
// The returned graph keeps each chunk's per-fragment path names; the caller
// runs merge_path_fragments once across all buckets. A bucket of one input is
// returned as-is (its graph is moved out, no stitching needed).
static unique_ptr<MutablePathDeletableHandleGraph>
stitch_reference_bucket(std::vector<MergeInput>& inputs, const Logger& logger, bool fusion) {

    using handlegraph::PathSense;
    using handlegraph::PathMetadata;
    using handlegraph::subrange_t;
    using handlegraph::offset_t;
    using handlegraph::nid_t;

    // (2) Sort by REFERENCE start offset.
    std::sort(inputs.begin(), inputs.end(),
              [](const MergeInput& a, const MergeInput& b) {
                  return a.ref_start_offset < b.ref_start_offset;
              });

    // (3) Validate that every chunk's paths share the appropriate boundary.
    for (size_t i = 0; i < inputs.size(); ++i) {
        auto& in = inputs[i];
        bool need_end_check = (i + 1 < inputs.size());
        bool need_start_check = (i > 0);
        in.graph->for_each_path_handle([&](const path_handle_t& p) {
            if (need_end_check) {
                step_handle_t last_step = in.graph->path_back(p);
                handle_t last_h = in.graph->get_handle_of_step(last_step);
                if (last_h != in.ref_end_handle) {
                    logger.error() << "Path \"" << in.graph->get_path_name(p)
                                   << "\" in \"" << in.filename
                                   << "\" does not end at the REFERENCE end node (id "
                                   << in.graph->get_id(in.ref_end_handle)
                                   << "); --merge/-m requires every path in a left/middle "
                                   << "chunk to terminate on the REFERENCE boundary." << endl;
                }
            }
            if (need_start_check) {
                step_handle_t first_step = in.graph->path_begin(p);
                handle_t first_h = in.graph->get_handle_of_step(first_step);
                if (first_h != in.ref_start_handle) {
                    logger.error() << "Path \"" << in.graph->get_path_name(p)
                                   << "\" in \"" << in.filename
                                   << "\" does not start at the REFERENCE start node (id "
                                   << in.graph->get_id(in.ref_start_handle)
                                   << "); --merge/-m requires every path in a right/middle "
                                   << "chunk to begin on the REFERENCE boundary." << endl;
                }
            }
        });
    }

    // (4) Pairwise stitch. inputs[0]'s graph becomes the destination accumulator.
    unique_ptr<MutablePathDeletableHandleGraph> dest = std::move(inputs[0].graph);
    int64_t max_node_id = dest->max_node_id();

    // Rightmost boundary handle, tracked in dest's ID space.
    handle_t dest_end_handle = inputs[0].ref_end_handle;
    offset_t dest_ref_end_offset = inputs[0].ref_end_offset;

    for (size_t i = 1; i < inputs.size(); ++i) {
        auto& in = inputs[i];
        unique_ptr<MutablePathDeletableHandleGraph> src = std::move(in.graph);

        // Capture src's REF endpoint identity BEFORE any mutation (divide /
        // increment may invalidate the original handles).
        nid_t  src_ref_start_id  = src->get_id(in.ref_start_handle);
        bool   src_ref_start_rev = src->get_is_reverse(in.ref_start_handle);
        nid_t  src_ref_end_id    = src->get_id(in.ref_end_handle);
        bool   src_ref_end_rev   = src->get_is_reverse(in.ref_end_handle);
        size_t src_start_len     = src->get_length(in.ref_start_handle);

        // Length of the left chunk's boundary node (only matters for fusion:
        // under -u every right-chunk path's first node grows by this many bp,
        // so the path metadata must shift back by this amount so the merged
        // subrange end still reflects the original genomic coordinate).
        size_t left_boundary_len = dest->get_length(dest_end_handle);

        // Overlap from REF offsets.
        offset_t overlap = dest_ref_end_offset - in.ref_start_offset;
        if (overlap < 0) {
            logger.error() << "REFERENCE offset gap between \"" << inputs[i-1].filename
                           << "\" (REF ends at offset " << dest_ref_end_offset
                           << ") and \"" << in.filename << "\" (REF starts at offset "
                           << in.ref_start_offset << "). --merge/-m does not insert gaps."
                           << endl;
        }

        // Trim the right chunk's start node by `overlap` bp.
        nid_t new_start_id_in_src = src_ref_start_id;
        if (overlap > 0) {
            if ((size_t)overlap > src_start_len) {
                logger.error() << "REFERENCE overlap of " << overlap << " bp between \""
                               << inputs[i-1].filename << "\" and \"" << in.filename
                               << "\" exceeds the right chunk's first node length ("
                               << src_start_len << " bp on node " << src_ref_start_id
                               << ")." << endl;
            }
            if ((size_t)overlap == src_start_len) {
                logger.error() << "REFERENCE overlap of " << overlap
                               << " bp exactly equals the right chunk's first node length; "
                               << "trim would yield a zero-length start node. Re-chunk so the "
                               << "first node of \"" << in.filename
                               << "\" extends past the overlap." << endl;
            }
            // Split the start node. divide_handle keeps paths consistent:
            // every path now visits piece_A then piece_B at its head.
            auto pieces = src->divide_handle(in.ref_start_handle,
                                             std::vector<size_t>{(size_t)overlap});
            handle_t piece_A = pieces[0];
            handle_t piece_B = pieces[1];
            new_start_id_in_src = src->get_id(piece_B);

            // Drop piece_A from every path's front. We materialize the path
            // list first to avoid iterator invalidation.
            std::vector<path_handle_t> all_paths;
            src->for_each_path_handle([&](const path_handle_t& p) {
                all_paths.push_back(p);
            });
            for (path_handle_t p : all_paths) {
                src->pop_front_step(p);
            }

            // piece_A has no remaining path steps and (because the boundary
            // node is the true left edge of the chunk) no inbound edges.
            // Its sole outbound edge piece_A->piece_B vanishes with the node.
            src->destroy_handle(piece_A);
        }

        // Reconcile each path's offset metadata with what its first node
        // actually represents in the destination accumulator:
        //  - Trimming removes `overlap` bp from the front, so every path
        //    "starts" `overlap` bp later in its own coordinate system.
        //  - Under --fusion, the right chunk's first node will gain
        //    `left_boundary_len` bp on its left (the left chunk's boundary
        //    sequence is fused in), so the path now "starts" that many bp
        //    earlier than just-after-trim.
        // The combined shift makes fragment_start + path_length come out to
        // the correct genomic end after fragment-merge.
        int64_t meta_shift = (int64_t)overlap - (fusion ? (int64_t)left_boundary_len : 0);
        if (meta_shift != 0) {
            struct PathSnapshot {
                PathSense sense;
                string sample;
                string locus;
                size_t haplotype;
                size_t phase_block;
                subrange_t subrange;
                bool is_circular;
                std::vector<handle_t> steps;
            };
            std::vector<PathSnapshot> snapshots;
            src->for_each_path_handle([&](const path_handle_t& sp) {
                PathSnapshot s;
                s.sense = src->get_sense(sp);
                s.sample = src->get_sample_name(sp);
                s.locus = src->get_locus_name(sp);
                s.haplotype = src->get_haplotype(sp);
                s.phase_block = src->get_phase_block(sp);
                s.subrange = src->get_subrange(sp);
                s.is_circular = src->get_is_circular(sp);
                for (handle_t h : src->scan_path(sp)) {
                    s.steps.push_back(h);
                }
                snapshots.push_back(std::move(s));
            });
            std::vector<path_handle_t> to_destroy;
            src->for_each_path_handle([&](const path_handle_t& sp) {
                to_destroy.push_back(sp);
            });
            for (path_handle_t sp : to_destroy) {
                src->destroy_path(sp);
            }
            for (auto& s : snapshots) {
                if (s.subrange != PathMetadata::NO_SUBRANGE) {
                    int64_t new_first = (int64_t)s.subrange.first + meta_shift;
                    if (new_first < 0) new_first = 0;
                    s.subrange.first = (offset_t)new_first;
                } else if (s.phase_block != PathMetadata::NO_PHASE_BLOCK) {
                    int64_t new_pb = (int64_t)s.phase_block + meta_shift;
                    if (new_pb < 0) new_pb = 0;
                    s.phase_block = (size_t)new_pb;
                } else {
                    int64_t new_first = meta_shift;
                    if (new_first < 0) new_first = 0;
                    s.subrange = subrange_t{(offset_t)new_first,
                                            PathMetadata::NO_END_POSITION};
                }
                path_handle_t new_p = src->create_path(
                    s.sense, s.sample, s.locus, s.haplotype,
                    s.phase_block, s.subrange, s.is_circular);
                for (handle_t h : s.steps) {
                    src->append_step(new_p, h);
                }
            }
        }

        // Renumber src out of dest's range.
        int64_t shift = 0;
        int64_t delta = max_node_id - (int64_t)src->min_node_id();
        if (delta >= 0) {
            src->increment_node_ids(delta + 1);
            shift = delta + 1;
        }
        max_node_id = std::max(max_node_id, (int64_t)src->max_node_id());

        // Copy nodes/edges, then each path under its original (per-fragment)
        // name. The fragment-merge pass at the end re-collapses by metadata.
        handlealgs::copy_handle_graph(src.get(), dest.get());
        src->for_each_path_handle([&](const path_handle_t& src_path) {
            string path_name = src->get_path_name(src_path);
            if (dest->has_path(path_name)) {
                logger.error() << "Path \"" << path_name
                               << "\" exists in multiple inputs to --merge/-m; "
                               << "fragments must carry distinct subrange-suffixed names."
                               << endl;
            }
            handlealgs::copy_path(src.get(), src_path, dest.get());
        });

        // Locate the right chunk's new start and REF end in dest's ID space.
        handle_t src_new_start_in_dest =
            dest->get_handle((nid_t)((int64_t)new_start_id_in_src + shift), src_ref_start_rev);
        handle_t src_ref_end_in_dest =
            dest->get_handle((nid_t)((int64_t)src_ref_end_id + shift), src_ref_end_rev);

        if (fusion) {
            // Glue dest_end_handle (left chunk's end) and src_new_start_in_dest
            // (right chunk's trimmed start) into a single node with the
            // concatenated forward-strand sequence. Both are visited only in
            // forward orientation by the paths that touch them (we asserted
            // forward boundaries above), so the concatenation order matches
            // what every path expects.
            string fused_seq = dest->get_sequence(dest_end_handle)
                             + dest->get_sequence(src_new_start_in_dest);
            handle_t fused = dest->create_handle(fused_seq);
            max_node_id = std::max(max_node_id, (int64_t)dest->get_id(fused));

            // Preserve every edge that wasn't the implicit seam.
            // - inbound to the left side of dest_end_handle -> inbound to fused
            // - outbound from the right side of src_new_start_in_dest ->
            //   outbound from fused
            std::vector<handle_t> preds;
            dest->follow_edges(dest_end_handle, true, [&](const handle_t& p) {
                preds.push_back(p);
            });
            std::vector<handle_t> succs;
            dest->follow_edges(src_new_start_in_dest, false, [&](const handle_t& n) {
                succs.push_back(n);
            });
            for (handle_t p : preds) {
                if (!dest->has_edge(p, fused)) dest->create_edge(p, fused);
            }
            for (handle_t n : succs) {
                if (!dest->has_edge(fused, n)) dest->create_edge(fused, n);
            }

            // Rewrite every step on the two originals to point at fused,
            // preserving the step's orientation defensively (we already
            // asserted forward but the code stays correct either way).
            auto replace_steps_on = [&](handle_t old_h) {
                std::vector<step_handle_t> steps;
                dest->for_each_step_on_handle(old_h, [&](const step_handle_t& s) {
                    steps.push_back(s);
                });
                for (const step_handle_t& s : steps) {
                    handle_t h = dest->get_handle_of_step(s);
                    handle_t replacement = dest->get_is_reverse(h)
                        ? dest->flip(fused) : fused;
                    step_handle_t next = dest->get_next_step(s);
                    dest->rewrite_segment(s, next, std::vector<handle_t>{replacement});
                }
            };
            replace_steps_on(dest_end_handle);
            replace_steps_on(src_new_start_in_dest);

            // Destroy originals (also removes any stray edges we didn't pick up).
            dest->destroy_handle(dest_end_handle);
            dest->destroy_handle(src_new_start_in_dest);

            // If the right chunk's REF is a single node, the new boundary IS
            // the fused node; otherwise it's the right chunk's REF end.
            dest_end_handle = (src_ref_end_id == new_start_id_in_src)
                ? fused
                : src_ref_end_in_dest;
        } else {
            if (!dest->has_edge(dest_end_handle, src_new_start_in_dest)) {
                dest->create_edge(dest_end_handle, src_new_start_in_dest);
            }
            dest_end_handle = src_ref_end_in_dest;
        }

        // The right chunk's REF end offset becomes the new dest end offset.
        // (Trimming only removed bases from the front, so the end is unchanged.)
        dest_ref_end_offset = in.ref_end_offset;

        logger.info() << "Stitched \"" << in.filename << "\": overlap=" << overlap
                      << " bp, " << (fusion ? "fused" : "edged") << " boundary."
                      << endl;
    }

    return dest;
}

// --merge / -m mode (optionally with --fusion / -u).
//
// Combines chunked-assembly graphs that share PathSense::REFERENCE paths
// (e.g. CHM13#0#chr22 / CHM13#0#chr22[17475777]) and overlap on them. The
// REFERENCE-path offsets pin each chunk to a coordinate; the overlap between
// adjacent chunks is computed from those offsets and trimmed off the right
// chunk's leading node before stitching.
//
// Inputs may describe more than one reference (e.g. some chunks of CHM13#0#chr21
// and some of CHM13#0#chr22). They are bucketed by REFERENCE identity (sample,
// locus, haplotype, circular); each bucket is stitched independently, then the
// per-bucket graphs are combined into one output (node IDs renumbered to avoid
// collisions, no edges added between distinct references). The final output
// carries one REFERENCE path per chromosome (e.g. both CHM13#0#chr21 and
// CHM13#0#chr22).
//
// Pipeline:
//   (1) Load all inputs and locate each graph's unique REFERENCE-sense path.
//       Require exactly one REFERENCE path per graph.
//   (2) Bucket inputs by REFERENCE identity.
//   (3) Stitch each bucket (see stitch_reference_bucket: sort, validate
//       boundaries, pairwise trim+connect/fuse), then copy the stitched bucket
//       graphs into a single destination, renumbering IDs as needed.
//   (4) Run the same fragment-merge pass as --connect-fragments to collapse
//       each reference's chunk paths into one coordinate-spanning path. Because
//       that pass buckets by (sense, sample, locus, haplotype, circular),
//       different chromosomes stay separate.
int merge_combine(int argc, char** argv, const Logger& logger, bool fusion) {

    using handlegraph::PathSense;
    using handlegraph::PathMetadata;
    using handlegraph::subrange_t;
    using handlegraph::offset_t;
    using handlegraph::nid_t;

    // (1) Load all inputs.
    std::vector<MergeInput> inputs;
    while (optind < argc) {
        MergeInput in;
        in.filename = get_input_file_name(optind, argc, argv);
        in.graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in.filename);
        inputs.push_back(std::move(in));
    }
    if (inputs.size() < 2) {
        logger.error() << "--merge/-m requires at least two input graphs." << endl;
    }

    // Find the unique REFERENCE path in each graph and pull its endpoint and
    // offset metadata.
    for (auto& in : inputs) {
        std::vector<path_handle_t> refs;
        in.graph->for_each_path_handle([&](const path_handle_t& p) {
            if (in.graph->get_sense(p) == PathSense::REFERENCE) {
                refs.push_back(p);
            }
        });
        if (refs.empty()) {
            logger.error() << "No REFERENCE-sense path found in \"" << in.filename
                           << "\"; --merge/-m requires exactly one per graph." << endl;
        }
        if (refs.size() > 1) {
            logger.error() << refs.size() << " REFERENCE-sense paths found in \""
                           << in.filename << "\"; --merge/-m requires exactly one." << endl;
        }
        in.ref_path = refs[0];
        in.ref_name = in.graph->get_path_name(in.ref_path);
        in.sample = in.graph->get_sample_name(in.ref_path);
        in.locus = in.graph->get_locus_name(in.ref_path);
        in.haplotype = in.graph->get_haplotype(in.ref_path);
        in.is_circular = in.graph->get_is_circular(in.ref_path);

        subrange_t sr = in.graph->get_subrange(in.ref_path);
        size_t pb = in.graph->get_phase_block(in.ref_path);
        if (sr != PathMetadata::NO_SUBRANGE) {
            in.ref_start_offset = sr.first;
        } else if (pb != PathMetadata::NO_PHASE_BLOCK) {
            in.ref_start_offset = (offset_t)pb;
        } else {
            in.ref_start_offset = 0;
        }
        size_t ref_len = compute_path_length(in.graph.get(), in.ref_path);
        in.ref_end_offset = in.ref_start_offset + (offset_t)ref_len;

        step_handle_t first_step = in.graph->path_begin(in.ref_path);
        step_handle_t last_step = in.graph->path_back(in.ref_path);
        in.ref_start_handle = in.graph->get_handle_of_step(first_step);
        in.ref_end_handle   = in.graph->get_handle_of_step(last_step);

        if (in.graph->get_is_reverse(in.ref_start_handle)
            || in.graph->get_is_reverse(in.ref_end_handle)) {
            logger.error() << "REFERENCE path \"" << in.ref_name << "\" in \""
                           << in.filename << "\" visits its boundary node in reverse "
                           << "orientation; --merge/-m only supports forward-strand boundaries."
                           << endl;
        }
    }

    // (2) Bucket inputs by REFERENCE identity (sample, locus, haplotype,
    //     circular). Chunks of the same chromosome land in one bucket; different
    //     chromosomes land in different buckets and are stitched separately.
    //     bucket_order preserves first-seen order for deterministic output.
    using BucketKey = std::tuple<string, string, size_t, bool>;
    std::map<BucketKey, std::vector<size_t>> buckets;
    std::vector<BucketKey> bucket_order;
    for (size_t i = 0; i < inputs.size(); ++i) {
        BucketKey key{inputs[i].sample, inputs[i].locus,
                      inputs[i].haplotype, inputs[i].is_circular};
        if (!buckets.count(key)) {
            bucket_order.push_back(key);
        }
        buckets[key].push_back(i);
    }

    // (3) Stitch each bucket, then combine the per-bucket graphs into a single
    //     destination accumulator. The first bucket's stitched graph is the
    //     destination; each subsequent one is renumbered out of the way and its
    //     nodes/edges/paths copied in. No edges are created between buckets, so
    //     distinct references remain disconnected components.
    unique_ptr<MutablePathDeletableHandleGraph> dest;
    int64_t max_node_id = 0;
    for (const auto& key : bucket_order) {
        std::vector<MergeInput> bucket;
        bucket.reserve(buckets[key].size());
        for (size_t i : buckets[key]) {
            bucket.push_back(std::move(inputs[i]));
        }

        unique_ptr<MutablePathDeletableHandleGraph> stitched =
            stitch_reference_bucket(bucket, logger, fusion);

        if (!dest) {
            dest = std::move(stitched);
            max_node_id = dest->max_node_id();
            continue;
        }

        int64_t delta = max_node_id - (int64_t)stitched->min_node_id();
        if (delta >= 0) {
            stitched->increment_node_ids(delta + 1);
        }
        max_node_id = std::max(max_node_id, (int64_t)stitched->max_node_id());

        handlealgs::copy_handle_graph(stitched.get(), dest.get());
        stitched->for_each_path_handle([&](const path_handle_t& src_path) {
            string path_name = stitched->get_path_name(src_path);
            if (dest->has_path(path_name)) {
                logger.error() << "Path \"" << path_name
                               << "\" exists in multiple --merge/-m references; inputs for "
                               << "different chromosomes must not share path names." << endl;
            }
            handlealgs::copy_path(stitched.get(), src_path, dest.get());
        });
    }

    // (4) Bucket and merge fragments into single coordinate-spanning paths
    //     (one per reference).
    merge_path_fragments(dest.get(), logger);

    // (5) Serialize.
    vg::io::save_handle_graph(dest.get(), cout);
    return 0;
}
