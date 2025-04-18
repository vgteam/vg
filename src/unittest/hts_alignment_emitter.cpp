/// \file hts_alignment_emitter.cpp
///  
/// unit tests for alignment emitter creation

#include <sys/stat.h>
#include <iostream>
#include "../handle.hpp"
#include "../utility.hpp"
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "catch.hpp"

#include "../hts_alignment_emitter.hpp"

namespace vg {
namespace unittest {

/// Make a graph and return a handle to a reference path
static path_handle_t populate_test_graph(bdsg::HashGraph& graph) {
    handle_t node_handle = graph.create_handle("GATTACA");
    path_handle_t ref_path = graph.create_path(PathSense::REFERENCE, "RefSamp", "chr1", 0, PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE, false);
    graph.append_step(ref_path, node_handle);
    return ref_path;
}

/// Make a sequence dictionary for the overlayed graph
static SequenceDictionary test_sequence_dictionary(const PathPositionHandleGraph* path_position_handle_graph, const bdsg::HashGraph& graph, const path_handle_t& ref_path) {
    // There's 1 reference path, and its graph length is the same as its base path length.
    path_handle_t overlayed_path = path_position_handle_graph->get_path_handle(graph.get_path_name(ref_path));
    size_t ref_length = path_position_handle_graph->get_path_length(overlayed_path);
    SequenceDictionary dict;
    dict.emplace_back();
    dict.back().path_handle = overlayed_path;
    dict.back().path_name = graph.get_path_name(ref_path);
    dict.back().path_length = ref_length;
    dict.back().base_path_name = dict.back().path_name;
    dict.back().base_path_length = dict.back().path_length;
    // TODO: Just use get_sdequence_dictionary instead?
    return dict;
}

/// Get the size of a file
static size_t get_file_size(const std::string& filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != 0) {
        throw std::runtime_error("Could not stat " + filename);
    }
    return buf.st_size;
}

/// Get a test unmapped read.
static Alignment get_unmapped_read() {
    // Write a read.
    Alignment unmapped;
    unmapped.set_sequence("NNN");
    unmapped.set_name("unmapped");
    // To prove this is unmapped in a linear reference and not just not surjected, it has to have *a* refpos, even if it's completely empty.
    unmapped.add_refpos();
    return unmapped;
}

TEST_CASE("Can create and use a BAM alignment emitter", "[giraffe][alignment_emitter][hts_alignment_emitter]") {
    // Set up the test fixtures
    bdsg::HashGraph graph;
    path_handle_t ph1 = populate_test_graph(graph);
    bdsg::ReferencePathOverlayHelper overlay_helper;
    PathPositionHandleGraph* overlay_graph = overlay_helper.apply(&graph);
    auto sequence_dictionary = test_sequence_dictionary(overlay_graph, graph, ph1);

    // Set up the temp file to write
    std::string out_dir = vg::temp_file::create_directory();
    std::string out_file = out_dir + "/test.bam";
    
    // Make the emitter
    unique_ptr<AlignmentEmitter> alignment_emitter = get_alignment_emitter(out_file, "BAM", sequence_dictionary, 1, overlay_graph, ALIGNMENT_EMITTER_FLAG_HTS_RAW);

    // We should get one
    REQUIRE(alignment_emitter);
    
    // Write a read.
    alignment_emitter->emit_single(get_unmapped_read());

    // Delete it to shut it down
    alignment_emitter.reset();

    // We should see data on disk
    REQUIRE(get_file_size(out_file) > 0);
}

TEST_CASE("Can create and use a CRAM alignment emitter", "[giraffe][alignment_emitter][hts_alignment_emitter]") {
    // Set up the test fixtures
    bdsg::HashGraph graph;
    path_handle_t ph1 = populate_test_graph(graph);
    bdsg::ReferencePathOverlayHelper overlay_helper;
    PathPositionHandleGraph* overlay_graph = overlay_helper.apply(&graph);
    auto sequence_dictionary = test_sequence_dictionary(overlay_graph, graph, ph1);

    // Set up the temp file to write
    std::string out_dir = vg::temp_file::create_directory();
    std::string out_file = out_dir + "/test.cram";
    
    // Make the emitter
    unique_ptr<AlignmentEmitter> alignment_emitter = get_alignment_emitter(out_file, "CRAM", sequence_dictionary, 1, overlay_graph, ALIGNMENT_EMITTER_FLAG_HTS_RAW);

    // We should get one
    REQUIRE(alignment_emitter);

    // Write a read.
    alignment_emitter->emit_single(get_unmapped_read());

    // Delete it to shut it down
    alignment_emitter.reset();

    // We should see data on disk
    REQUIRE(get_file_size(out_file) > 0);
}



}
}
    
