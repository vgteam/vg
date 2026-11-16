/** \file
 *
 * Unit tests for recombinator.cpp: clipping gref fragments to a sampled graph,
 * selecting haplotypes in a subchain, and threading the generated haplotypes.
 */

#include "../recombinator.hpp"

#include "catch.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// gbwt::vector_type may store nodes in a narrower type than gbwt::node_type.
typedef gbwt::vector_type::value_type node_t;
node_t fwd(nid_t id) { return gbwt::Node::encode(id, false); }
node_t rev(nid_t id) { return gbwt::Node::encode(id, true); }

// Node length is the node id, so offsets are easy to predict.
size_t id_as_length(gbwt::node_type node) { return gbwt::Node::id(node); }

void check_piece(const PathPiece& piece, size_t start, size_t nodes, size_t bp_offset, size_t bp_length) {
    REQUIRE(piece.start == start);
    REQUIRE(piece.nodes == nodes);
    REQUIRE(piece.bp_offset == bp_offset);
    REQUIRE(piece.bp_length == bp_length);
}

} // anonymous namespace

//------------------------------------------------------------------------------

TEST_CASE("Paths are clipped to the nodes and edges of a GBWT index", "[gref][recombinator]") {
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    // Paths 1-2-3, 3-4, 5-6 and 8. Node 7 is inside the alphabet range but unused,
    // and nodes 4 and 5 are both present with no edge between them.
    gbwt::GBWTBuilder builder(gbwt::bit_length(rev(8)), 1024);
    std::vector<gbwt::vector_type> paths {
        { fwd(1), fwd(2), fwd(3) }, { fwd(3), fwd(4) }, { fwd(5), fwd(6) }, { fwd(8) }
    };
    for (const gbwt::vector_type& path : paths) {
        builder.insert(path, true);
    }
    builder.finish();
    const gbwt::DynamicGBWT& index = builder.index;
    REQUIRE(index.contains(fwd(7)));

    SECTION("a path that is entirely in the index is one piece") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(3) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("the reverse orientation is in the index too") {
        gbwt::vector_type path { rev(3), rev(2), rev(1) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 3, 0, 6);
    }

    SECTION("a missing node splits the path and still counts toward the offset") {
        gbwt::vector_type path { fwd(1), fwd(2), fwd(7), fwd(3), fwd(4) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 3);
        check_piece(pieces[1], 3, 2, 10, 7);
    }

    SECTION("present nodes without an edge between them are split") {
        gbwt::vector_type path { fwd(3), fwd(4), fwd(5), fwd(6) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 2);
        check_piece(pieces[0], 0, 2, 0, 7);
        check_piece(pieces[1], 2, 2, 7, 11);
    }

    SECTION("a node outside the alphabet range is missing") {
        gbwt::vector_type path { fwd(2), fwd(100) };
        auto pieces = clip_path_to_index(path, index, id_as_length);
        REQUIRE(pieces.size() == 1);
        check_piece(pieces[0], 0, 1, 0, 2);
    }

    SECTION("nothing survives in an empty index") {
        gbwt::DynamicGBWT empty;
        gbwt::vector_type path { fwd(1), fwd(2) };
        REQUIRE(clip_path_to_index(path, empty, id_as_length).empty());
    }
}

//------------------------------------------------------------------------------

namespace {

// A GBZ graph with one subchain from node 1 to node 6 and four haplotypes:
// sequence 0: 1 2 4 6
// sequence 1: 1 3 4 6 (node 3 is 2000 bp)
// sequence 2: 1 2 5 6 (node 5 has 6 N bases)
// sequence 3: 1 2 4 6
struct SelectionFixture {
    gbwtgraph::GBZ gbz;
    gbwt::FragmentMap fragment_map;
    Haplotypes::Subchain subchain;
    hash_map<Haplotypes::Subchain::kmer_type, size_t> kmer_counts;
    double coverage = 30.0;

    SelectionFixture() {
        std::vector<std::vector<nid_t>> paths { { 1, 2, 4, 6 }, { 1, 3, 4, 6 }, { 1, 2, 5, 6 }, { 1, 2, 4, 6 } };
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
        gbwt::GBWTBuilder builder(gbwt::bit_length(gbwt::Node::encode(6, true)), 100);
        for (const auto& path : paths) {
            gbwt::vector_type encoded;
            for (nid_t id : path) {
                encoded.push_back(gbwt::Node::encode(id, false));
            }
            builder.insert(encoded, true);
        }
        builder.finish();
        builder.index.addMetadata();
        builder.index.metadata.setSamples({ "s0", "s1", "s2", "s3" });
        builder.index.metadata.setHaplotypes(paths.size());
        builder.index.metadata.setContigs({ "chr" });
        for (size_t i = 0; i < paths.size(); i++) {
            builder.index.metadata.addPath(i, 0, 0, 0);
        }
        gbwt::GBWT index(builder.index);

        gbwtgraph::NaiveGraph source;
        source.create_node(1, "GATTACAGATTACA");
        source.create_node(2, "ACGTA");
        std::string long_node;
        for (size_t i = 0; i < 500; i++) {
            long_node += "ACGT";
        }
        source.create_node(3, long_node);
        source.create_node(4, "CCCGG");
        source.create_node(5, "ACNNNNNNGT");
        source.create_node(6, "TTGACCA");
        this->gbz = gbwtgraph::GBZ(index, source);
        this->fragment_map = gbwt::FragmentMap(this->gbz.index.metadata, false);

        // Kmers 0 to 3 are in sequences 0 and 3, kmers 4 and 5 in sequence 1, and
        // kmer 6 in sequence 2. All of them are present in the sample.
        this->subchain.type = Haplotypes::Subchain::normal;
        this->subchain.start = gbwt::Node::encode(1, false);
        this->subchain.end = gbwt::Node::encode(6, false);
        this->subchain.placeholder = 0;
        for (size_t i = 0; i < paths.size(); i++) {
            // Sequence `i` starts at offset `i` of node 1.
            this->subchain.sequences.push_back({ gbwt::Path::encode(i, false), i });
        }
        this->set_rows({ "1111000", "0000110", "0000001", "1111000" });
        for (size_t j = 0; j < this->subchain.num_kmers(); j++) {
            this->kmer_counts[j + 1] = 30;
        }
    }

    // Sets the kmers of the subchain: kmer `j` is in sequence `i` if `rows[i][j]` is `1`.
    void set_rows(const std::vector<std::string>& rows) {
        size_t num_kmers = rows.front().size();
        sdsl::bit_vector matrix(rows.size() * num_kmers, 0);
        for (size_t i = 0; i < rows.size(); i++) {
            for (size_t j = 0; j < num_kmers; j++) {
                matrix[i * num_kmers + j] = (rows[i][j] == '1');
            }
        }
        this->subchain.kmers.clear();
        for (size_t j = 0; j < num_kmers; j++) {
            this->subchain.kmers.push_back(j + 1);
        }
        this->subchain.kmer_counts = sdsl::int_vector<0>(num_kmers, 2, 64);
        this->subchain.kmers_present = KmerPresenceMatrix(rows.size(), num_kmers, std::move(matrix));
    }

    // Returns the sequence offsets of the selected haplotypes.
    std::vector<size_t> select(const Recombinator::Parameters& parameters, const Recombinator::SelectionContext& context) const {
        Recombinator::Statistics statistics;
        auto selected = select_haplotypes(
            this->gbz, this->subchain, this->kmer_counts, this->coverage, 29, &statistics, nullptr, parameters, context
        );
        std::vector<size_t> result;
        for (auto& haplotype : selected) {
            result.push_back(haplotype.first);
        }
        return result;
    }
};

} // anonymous namespace

TEST_CASE("Haplotype selection in a subchain", "[haplotypes][recombinator]") {
    SelectionFixture fixture;
    Recombinator::Parameters parameters;
    parameters.num_haplotypes = 4;
    Recombinator::SelectionContext context;
    context.fragment_map = &fixture.fragment_map;

    SECTION("the fixture is what the test expects") {
        for (size_t i = 0; i < 4; i++) {
            gbwt::edge_type pos(fixture.subchain.start, fixture.subchain.sequences[i].second);
            gbwt::edge_type next = fixture.gbz.index.LF(pos);
            REQUIRE(gbwt::Node::id(next.first) == (i == 1 ? 3 : 2));
        }
    }

    SECTION("without limits") {
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 1, 2 }));
    }

    SECTION("new-bp limit") {
        parameters.max_new_bp = 1000;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 2, 0 }));
    }

    SECTION("new-bp limit with an escape") {
        parameters.max_new_bp = 1000;
        parameters.bp_escape = 1000.0; // 2 new sample kmers for 2000 bp.
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 1, 2 }));
        parameters.bp_escape = 999.0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 2, 0 }));
    }

    SECTION("the escape counts only sample kmers") {
        fixture.kmer_counts[5] = 0; // Kmer 4 is absent from the sample.
        parameters.max_new_bp = 1000;
        parameters.bp_escape = 1000.0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 2, 0 }));
        parameters.bp_escape = 2000.0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 2, 1 }));
    }

    SECTION("the escape counts only kmers that the selected haplotypes lack") {
        // Sequence 3 also has kmer 4, and it becomes the first haplotype.
        fixture.set_rows({ "1111000", "0000110", "0000001", "1111100" });
        parameters.max_new_bp = 1000;
        parameters.bp_escape = 1000.0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 3, 0, 2, 3 }));
        parameters.bp_escape = 2000.0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 3, 0, 1, 2 }));
    }

    SECTION("the first haplotype is exempt from the new-bp limit") {
        std::vector<gbwt::size_type> sources { fixture.subchain.sequences[1].first };
        context.sources = &sources;
        parameters.continuity_bonus = 10.0;
        parameters.max_new_bp = 1000;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 1, 0, 3, 2 }));
    }

    SECTION("with no candidate within the N-base limit, the best one is the first haplotype") {
        parameters.banned_samples = { "s0", "s1", "s3" };
        parameters.max_n_bp = 0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 2, 2, 2, 2 }));
    }

    SECTION("N-base limit") {
        parameters.max_n_bp = 0;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 1, 0 }));
        parameters.max_n_bp = 6;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 3, 1, 2 }));
    }

    SECTION("continuity bonus") {
        std::vector<gbwt::size_type> sources { fixture.subchain.sequences[2].first, gbwt::invalid_sequence() };
        context.sources = &sources;
        parameters.continuity_bonus = 10.0;
        REQUIRE(fixture.select(parameters, context).front() == 2);
    }

    SECTION("the model pair comes first and is exempt from the limits") {
        parameters.max_new_bp = 1000;
        parameters.max_n_bp = 0;
        std::pair<size_t, size_t> pair(2, 1);
        context.pair = &pair;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 2, 1, 0, 3 }));
    }

    SECTION("a haplotype may pair with itself") {
        std::pair<size_t, size_t> pair(0, 0);
        context.pair = &pair;
        REQUIRE(fixture.select(parameters, context) == std::vector<size_t>({ 0, 0, 3, 1 }));
    }

    SECTION("the model pair must be selectable") {
        parameters.banned_samples.insert("s1");
        std::pair<size_t, size_t> pair(0, 1);
        context.pair = &pair;
        REQUIRE_THROWS_AS(fixture.select(parameters, context), std::runtime_error);
    }
}

TEST_CASE("Threading by similarity", "[haplotypes][recombinator]") {
    // Sequence `i` is GBWT sequence `2 * i`.
    std::vector<std::string> rows { "1100", "1110", "0011", "0001" };
    Haplotypes::Subchain subchain;
    subchain.type = Haplotypes::Subchain::normal;
    sdsl::bit_vector matrix(rows.size() * 4, 0);
    for (size_t i = 0; i < rows.size(); i++) {
        subchain.sequences.push_back({ 2 * i, 0 });
        for (size_t j = 0; j < 4; j++) {
            matrix[i * 4 + j] = (rows[i][j] == '1');
        }
    }
    subchain.kmers_present = KmerPresenceMatrix(rows.size(), 4, std::move(matrix));
    std::vector<std::pair<size_t, double>> selected { { 1, 0.0 }, { 3, 0.0 } };

    SECTION("an unmatched haplotype takes the closest unused selection") {
        // Sources: sequence 2 (closest to 3) and sequence 0 (closest to 1).
        std::vector<gbwt::size_type> sources { 4, 0 };
        std::vector<size_t> haplotype_to_selected(2, 2);
        sdsl::bit_vector selected_in_use(2, 0);
        thread_by_similarity(subchain, selected, sources, haplotype_to_selected, selected_in_use);
        REQUIRE(haplotype_to_selected == std::vector<size_t>({ 1, 0 }));
        REQUIRE(selected_in_use[0] == 1);
        REQUIRE(selected_in_use[1] == 1);
    }

    SECTION("matched haplotypes and used selections are kept") {
        std::vector<gbwt::size_type> sources { 0, 4 };
        std::vector<size_t> haplotype_to_selected { 2, 0 };
        sdsl::bit_vector selected_in_use(2, 0);
        selected_in_use[0] = 1;
        thread_by_similarity(subchain, selected, sources, haplotype_to_selected, selected_in_use);
        REQUIRE(haplotype_to_selected == std::vector<size_t>({ 1, 0 }));
    }

    SECTION("repeated selections are used from the lowest index") {
        std::vector<std::pair<size_t, double>> repeated { { 1, 0.0 }, { 1, 0.0 }, { 3, 0.0 } };
        std::vector<gbwt::size_type> sources { 0, 0, gbwt::invalid_sequence() };
        std::vector<size_t> haplotype_to_selected(3, 3);
        sdsl::bit_vector selected_in_use(3, 0);
        thread_by_similarity(subchain, repeated, sources, haplotype_to_selected, selected_in_use);
        REQUIRE(haplotype_to_selected == std::vector<size_t>({ 0, 1, 3 }));
    }

    SECTION("a haplotype without a source in the subchain stays unmatched") {
        std::vector<gbwt::size_type> sources { 100, gbwt::invalid_sequence() };
        std::vector<size_t> haplotype_to_selected(2, 2);
        sdsl::bit_vector selected_in_use(2, 0);
        thread_by_similarity(subchain, selected, sources, haplotype_to_selected, selected_in_use);
        REQUIRE(haplotype_to_selected == std::vector<size_t>({ 2, 2 }));
        REQUIRE(selected_in_use[0] == 0);
    }
}

TEST_CASE("Model pair candidates in a subchain", "[haplotypes][recombinator]") {
    SelectionFixture fixture;
    Recombinator::Parameters parameters;
    pair_model::CountModel model;
    model.mu = 15.0; model.eps = 0.5; model.r1 = 20.0; model.r0 = 0.5;
    std::vector<double> table = pair_model::count_table(model, std::vector<double>(pair_model::CountModel::MAX_COUNT + 1, 1.0));
    auto candidates = [&]() -> pair_model::Candidates {
        return model_candidates(fixture.gbz, fixture.subchain, fixture.kmer_counts, table, fixture.fragment_map, 29, parameters);
    };
    // Returns the group that contains the sequence, or the number of groups.
    auto group_of = [](const pair_model::Candidates& result, size_t sequence) -> size_t {
        for (size_t g = 0; g < result.size(); g++) {
            if (std::find(result.members[g].begin(), result.members[g].end(), sequence) != result.members[g].end()) {
                return g;
            }
        }
        return result.size();
    };

    SECTION("groups, haplotypes, and the length of the subchain") {
        pair_model::Candidates result = candidates();
        REQUIRE(result.size() == 3);
        size_t g = group_of(result, 0);
        REQUIRE(result.members[g] == std::vector<size_t>({ 0, 3 }));
        // The haplotype of a sequence is (sample << 32) | phase.
        REQUIRE(result.haplotypes[g] == std::vector<std::uint64_t>({ 0, std::uint64_t(3) << 32 }));
        REQUIRE(group_of(result, 2) < result.size());
        REQUIRE(result.scores.size() == result.size() * result.size());
        // Sequence lengths are 31, 2026, 36, and 31 bp.
        REQUIRE(result.length == Approx(33.5));
    }

    SECTION("a sequence over the N-base limit is not a candidate") {
        parameters.max_n_bp = 5;
        pair_model::Candidates result = candidates();
        REQUIRE(result.size() == 2);
        REQUIRE(group_of(result, 2) == result.size());
        parameters.max_n_bp = 6;
        REQUIRE(group_of(candidates(), 2) < candidates().size());
    }

    SECTION("unless every sequence that can be selected is over the limit") {
        parameters.max_n_bp = 0;
        parameters.banned_samples = { "s0", "s1", "s3" };
        pair_model::Candidates result = candidates();
        REQUIRE(result.size() == 1);
        REQUIRE(result.members[0] == std::vector<size_t>({ 2 }));
    }

    SECTION("if no sequence can be selected, every sequence is a candidate") {
        parameters.banned_samples = { "s0", "s1", "s2", "s3" };
        pair_model::Candidates result = candidates();
        REQUIRE(result.size() == 3);
        for (size_t i = 0; i < fixture.subchain.num_sequences(); i++) {
            REQUIRE(group_of(result, i) < result.size());
        }
    }
}

TEST_CASE("Validating the sampling parameters", "[haplotypes][recombinator]") {
    Recombinator::Parameters parameters;
    REQUIRE_NOTHROW(validate_recombinator_parameters(parameters));

    SECTION("the bp escape needs a new-bp limit") {
        parameters.bp_escape = 300.0;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(parameters), std::runtime_error);
        parameters.max_new_bp = 1000;
        REQUIRE_NOTHROW(validate_recombinator_parameters(parameters));
    }

    SECTION("negative values") {
        Recombinator::Parameters bonus = parameters;
        bonus.continuity_bonus = -1.0;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(bonus), std::runtime_error);
        Recombinator::Parameters escape = parameters;
        escape.max_new_bp = 1000;
        escape.bp_escape = -1.0;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(escape), std::runtime_error);
    }

    SECTION("the model pair") {
        parameters.model_pair = true;
        parameters.num_haplotypes = 32;
        REQUIRE_NOTHROW(validate_recombinator_parameters(parameters));
        Recombinator::Parameters one = parameters;
        one.num_haplotypes = 1;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(one), std::runtime_error);
        Recombinator::Parameters diploid = parameters;
        diploid.diploid_sampling = true;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(diploid), std::runtime_error);
        Recombinator::Parameters haploid = parameters;
        haploid.haploid_scoring = true;
        REQUIRE_THROWS_AS(validate_recombinator_parameters(haploid), std::runtime_error);
    }

    SECTION("the model preset") {
        auto preset = Recombinator::Parameters::parse_preset("model");
        REQUIRE(preset);
        REQUIRE(*preset == Recombinator::Parameters::preset_model);
        Recombinator::Parameters model(*preset);
        REQUIRE(model.num_haplotypes == Recombinator::MODEL_NUM_HAPLOTYPES);
        REQUIRE(model.include_reference);
        REQUIRE(model.continuity_bonus == Recombinator::MODEL_CONTINUITY_BONUS);
        REQUIRE(model.max_new_bp == Recombinator::MODEL_MAX_NEW_BP);
        REQUIRE(model.bp_escape == Recombinator::MODEL_BP_ESCAPE);
        REQUIRE(model.max_n_bp == Recombinator::MODEL_MAX_N_BP);
        REQUIRE(model.similar_threading);
        REQUIRE(model.model_pair);
        REQUIRE_FALSE(model.diploid_sampling);
        REQUIRE_FALSE(model.haploid_scoring);
        REQUIRE_NOTHROW(validate_recombinator_parameters(model));
    }
}

//------------------------------------------------------------------------------

} // namespace unittest

} // namespace vg
