/// \file flat_records.cpp
///
/// Unit tests for writing a site as flat records: the base-level differences of its strands,
/// with the sites nested in it spliced in, from its resolved reference.

#include <set>

#include <bdsg/overlays/path_position_overlays.hpp>

#include "../../src/unittest/catch.hpp"
#include "../block_records.hpp"
#include "support/staged_nested_sites.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

/// The flat records of the top site of a `StagedNestedSites`, and what they were cut from.
struct TopSiteFlatRecords {
    /// The strands aligned to the reference, as a run of the one site, or nothing if they could
    /// not be.
    optional<FlatSite> flat;
    /// The records, in the order they were written.
    vector<vcflib::Variant> records;
    /// The number of lines the writer says it wrote.
    int written = 0;
};

/// Write the top site of `fixture` as flat records with `writer`, from a site record with the GT
/// `site_gt` whose first reference base is at position 101.
static TopSiteFlatRecords flatten_top_site(StagedNestedSites& fixture,
                                           const BlockRecordWriter& writer,
                                           const string& site_gt) {
    const StagedSite& top = fixture.site(fixture.top_key);
    TopSiteFlatRecords out;
    // The strands are spelled along the site's reference walk, as those of a run are.
    vector<string> strands;
    vector<int> travs;
    for (int strand = 0; strand < 2; ++strand) {
        strands.push_back(*fixture.splicer.sequence_along(
            top.alleles.walks[top.alleles.ref_allele], {&top}, strand, false));
        REQUIRE(strands.back() == *fixture.splicer.strand_sequence(top, strand));
        travs.push_back(*fixture.splicer.strand_traversal(top, strand));
    }
    FlatRefusal refused{};
    const shared_ptr<FlatRun> run =
        writer.align_flat(*fixture.splicer.reference_sequence(top), std::move(strands), refused);
    if (!run) {
        REQUIRE(refused == FlatRefusal::AlignerLimit);
        return out;
    }
    run->id = top.id;
    out.flat = FlatSite{run, 0, travs, 0, run->cuts.size()};
    SiteRecord record{};
    record.variant.sequenceName = "ref";
    record.variant.position = 101;
    record.variant.id = top.id;
    record.variant.ref = run->reference;
    record.variant.alleles = {run->reference};
    record.variant.format = {"GT", "PS"};
    record.variant.samples["sample"]["GT"] = {site_gt};
    record.variant.samples["sample"]["PS"] = {"7"};
    record.unflattened_position = 101;
    record.alleles_merged = false;
    for (const int trav : travs) {
        record.trav_to_allele[trav] = 0;
    }
    bdsg::PositionOverlay positioned(&fixture.nested.graph);
    out.written = writer.write_flat(positioned, *out.flat, "sample", nullptr, record,
                                    GLLayout::IMajor, [&](vcflib::Variant& line, size_t) {
                                        out.records.push_back(line);
                                        return true;
                                    });
    return out;
}

TEST_CASE("Flat records put into the resolved reference spell each strand", "[flat_records]") {
    StagedNestedSites fixture;
    BlockRecordWriter writer;

    SECTION("Differences of the two strands that share reference bases are a record each") {
        // Strand 0 leaves out node 4 at the inner site. Strand 1 takes node 6 at the middle
        // site, in place of nodes 3 to 5.
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.phase(fixture.inner_key, 1, -1, 1, 0);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "0|0");
        REQUIRE(top.flat);
        REQUIRE(top.written == 2);
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat->run, top.records, 101, "sample"));
        // Strand 1 has T for eight bases, and strand 0 lacks four of them. A strand has allele
        // 0 in the record of the other strand's difference.
        const vcflib::Variant& line = top.records[0];
        REQUIRE(line.position == 105);
        REQUIRE(line.ref == "GCTGAGCA");
        REQUIRE(line.alt == vector<string>{"T"});
        REQUIRE(line.samples.at("sample").at("GT") == vector<string>{"0|1"});
        REQUIRE(line.samples.at("sample").at("PS") == vector<string>{"7"});
        REQUIRE(line.id == ">1>8_0");
        REQUIRE(line.info.at("SB") == vector<string>{"0", "2"});
        REQUIRE(line.infoFlags.at("FL"));
        // No walk spells a flat allele.
        REQUIRE(line.info.at("AT") == vector<string>{".", "."});
        REQUIRE(top.records[1].position == 105);
        REQUIRE(top.records[1].ref == "GCTGA");
        REQUIRE(top.records[1].alt == vector<string>{"G"});
        REQUIRE(top.records[1].samples.at("sample").at("GT") == vector<string>{"1|0"});
        REQUIRE(top.records[1].id == ">1>8_1");

        // The records stand in for those of every site staged under the top site.
        set<size_t> under;
        fixture.splicer.for_each_under(fixture.top_key, [&](const StagedSite& nested) {
            under.insert(nested.record_key);
        });
        REQUIRE(under == set<size_t>{fixture.middle_key, fixture.inner_key, fixture.edit_key});
    }
    SECTION("A deletion on both strands is homozygous, with the base before it") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 0);
        fixture.phase(fixture.inner_key, 1, 1);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "0|0");
        REQUIRE(top.written == 1);
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat->run, top.records, 101, "sample"));
        REQUIRE(top.records[0].position == 105);
        REQUIRE(top.records[0].ref == "GCTGA");
        REQUIRE(top.records[0].alt == vector<string>{"G"});
        REQUIRE(top.records[0].samples.at("sample").at("GT") == vector<string>{"1|1"});
    }
    SECTION("A change inside a node on one strand is a heterozygous record") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 0);
        fixture.phase(fixture.inner_key, 0, 0);
        fixture.phase(fixture.edit_key, 0, 1);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "0|0");
        REQUIRE(top.written == 1);
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat->run, top.records, 101, "sample"));
        REQUIRE(top.records[0].position == 107);
        REQUIRE(top.records[0].ref == "T");
        REQUIRE(top.records[0].alt == vector<string>{"A"});
        REQUIRE(top.records[0].samples.at("sample").at("GT") == vector<string>{"0|1"});

        // The check fails for records that spell something else.
        vector<vcflib::Variant> wrong = top.records;
        wrong[0].samples["sample"]["GT"] = {"1|0"};
        REQUIRE(!BlockRecordWriter::flat_records_spell(*top.flat->run, wrong, 101, "sample"));
    }
    SECTION("Separate differences are separate records, unphased where the site is") {
        // Strand 0 deletes the middle site's nodes, from node 1 to node 8. Strand 1 has the
        // change inside node 4.
        fixture.phase(fixture.top_key, 1, 0);
        fixture.phase(fixture.middle_key, 0, -1, 1, 1);
        fixture.phase(fixture.inner_key, 0, -1, 1, 1);
        fixture.phase(fixture.edit_key, 1, -1, 1, 1);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "1/0");
        REQUIRE(top.flat);
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat->run, top.records, 101, "sample"));
        // The deletion covers the changed base, and each is a record of its strand.
        REQUIRE(top.written == 2);
        REQUIRE(top.records[0].samples.at("sample").at("GT") == vector<string>{"1/0"});
        REQUIRE(top.records[1].samples.at("sample").at("GT") == vector<string>{"0/1"});
        REQUIRE(!top.records[0].samples.at("sample").contains("PS"));
    }
    SECTION("Strands that spell the reference write nothing") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 0);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "0|0");
        REQUIRE(top.flat);
        REQUIRE(top.written == 0);
        REQUIRE(top.records.empty());
    }
    SECTION("A site the aligner cannot align within its limit is not flattened") {
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.index();
        StringAlignerParams tiny;
        tiny.max_cells = 4;
        writer.set_flat_aligner(tiny);
        REQUIRE(!flatten_top_site(fixture, writer, "0|0").flat);
    }
}

TEST_CASE("Flat records of a run of two sites are cut from one alignment of each strand",
          "[flat_records]") {
    BlockRecordWriter writer;
    // Strand 0 has a mismatch and an insertion in a repeat. Strand 1 has the same mismatch and a
    // deletion further on.
    const string reference = "GATTACAGGCACACACATTGCCGGAATT";
    const string first = "GATTCCAGGCACACACACATTGCCGGAATT";
    const string second = "GATTCCAGGCACACACATTGCAATT";
    FlatRefusal refused{};
    const shared_ptr<FlatRun> run = writer.align_flat(reference, {first, second}, refused);
    REQUIRE(run);
    REQUIRE(run->cuts.size() == 3);
    bdsg::HashGraph graph;
    bdsg::PositionOverlay positioned(&graph);
    // The record of a site whose first base is at `position`, with the GT `gt`.
    auto site_record = [&](const string& id, int64_t position, const string& gt) {
        SiteRecord record{};
        record.variant.sequenceName = "ref";
        record.variant.position = position;
        record.variant.id = id;
        record.variant.ref = "N";
        record.variant.alleles = {"N"};
        record.variant.format = {"GT"};
        record.variant.samples["sample"]["GT"] = {gt};
        record.unflattened_position = position;
        record.alleles_merged = false;
        record.trav_to_allele = {{0, 0}, {1, 1}, {2, 2}};
        return record;
    };
    vector<vcflib::Variant> records;
    auto keep = [&](vcflib::Variant& line, size_t) {
        records.push_back(line);
        return true;
    };

    SECTION("One site writes every record of its run, under its own ID") {
        run->id = ">1>2";
        const FlatSite flat{run, 0, {1, 2}, 0, 3};
        REQUIRE(writer.write_flat(positioned, flat, "sample", nullptr,
                                  site_record(">1>2", 1, "1|2"), GLLayout::IMajor, keep)
                == 3);
        REQUIRE(BlockRecordWriter::flat_records_spell(*run, records, 1, "sample"));
        // The mismatch, on both strands.
        REQUIRE(records[0].position == 5);
        REQUIRE(records[0].ref == "A");
        REQUIRE(records[0].alt == vector<string>{"C"});
        REQUIRE(records[0].samples.at("sample").at("GT") == vector<string>{"1|1"});
        // The insertion, at the start of the repeat, after the base before it.
        REQUIRE(records[1].position == 9);
        REQUIRE(records[1].ref == "G");
        REQUIRE(records[1].alt == vector<string>{"GCA"});
        REQUIRE(records[1].samples.at("sample").at("GT") == vector<string>{"1|0"});
        // The deletion.
        REQUIRE(records[2].ref.size() == 4);
        REQUIRE(records[2].alt[0].size() == 1);
        REQUIRE(records[2].samples.at("sample").at("GT") == vector<string>{"0|1"});
        REQUIRE(records[2].id == ">1>2_2");
        REQUIRE(records[2].info.at("SB") == vector<string>{"2", "3"});
    }
    SECTION("Each site of a run writes the records that start in it, numbered along the run") {
        // The second site's first base is the reference's base 13, at position 13, and each
        // record is placed from the first base of the site that writes it.
        run->id = ">1>3";
        const FlatSite before{run, 0, {1, 2}, 0, 2};
        const FlatSite after{run, 12, {0, 1}, 2, 1};
        REQUIRE(writer.write_flat(positioned, after, "sample", nullptr,
                                  site_record(">2>3", 13, "0|1"), GLLayout::IMajor, keep)
                == 1);
        REQUIRE(writer.write_flat(positioned, before, "sample", nullptr,
                                  site_record(">1>2", 1, "1|2"), GLLayout::IMajor, keep)
                == 2);
        std::ranges::sort(records, {}, &vcflib::Variant::position);
        REQUIRE(BlockRecordWriter::flat_records_spell(*run, records, 1, "sample"));
        REQUIRE(records[0].position == 5);
        REQUIRE(records[1].position == 9);
        REQUIRE(records[2].ref.size() == 4);
        REQUIRE(records[2].samples.at("sample").at("GT") == vector<string>{"0|1"});
        for (size_t i = 0; i < records.size(); ++i) {
            REQUIRE(records[i].id == ">1>3_" + std::to_string(i));
            REQUIRE(records[i].info.at("SB") == vector<string>{std::to_string(i), "3"});
        }
        // One site's records alone do not spell the strands.
        records.pop_back();
        REQUIRE(!BlockRecordWriter::flat_records_spell(*run, records, 1, "sample"));
    }
    SECTION("Strands are aligned at the costs of an alignment of two assemblies") {
        // The strand has GGG five bases before where the reference has it. Five of the eight
        // columns between would be mismatches, which cost 200, and the two gaps cost 198.
        const string moved_reference = "TTCATCAACGGAGGGTCATTCA";
        const string moved = "TTCATCAGGGACGGATCATTCA";
        const shared_ptr<FlatRun> gaps = writer.align_flat(moved_reference, {moved}, refused);
        REQUIRE(gaps);
        REQUIRE(gaps->strands[0].differences.size() == 2);
        REQUIRE(gaps->strands[0].differences[0].alt == "GGG");
        REQUIRE(gaps->strands[0].differences[0].ref_begin == 7);
        REQUIRE(gaps->strands[0].differences[0].ref_end == 7);
        REQUIRE(gaps->strands[0].differences[1].alt.empty());
        // At 4 for a mismatch and 6 + 2L for a gap, the mismatches cost less.
        writer.set_flat_aligner(StringAlignerParams());
        const shared_ptr<FlatRun> mismatches =
            writer.align_flat(moved_reference, {moved}, refused);
        REQUIRE(mismatches);
        REQUIRE(mismatches->strands[0].differences.size() == 3);
    }
    SECTION("A gap at the first base of a run has no base before it") {
        REQUIRE(!writer.align_flat("ACGT", {"CGT", "ACGT"}, refused));
        REQUIRE(refused == FlatRefusal::NoAnchorBase);
    }
}

}
}
