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
    /// The strands aligned to the resolved reference, or nothing if they could not be.
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
    vector<pair<int, string>> strands;
    for (int strand = 0; strand < 2; ++strand) {
        strands.emplace_back(*fixture.splicer.strand_traversal(top, strand),
                             *fixture.splicer.strand_sequence(top, strand));
    }
    out.flat = writer.align_flat(*fixture.splicer.reference_sequence(top), std::move(strands));
    if (!out.flat) {
        return out;
    }
    SiteRecord record{};
    record.variant.sequenceName = "ref";
    record.variant.position = 101;
    record.variant.id = top.id;
    record.variant.ref = out.flat->reference;
    record.variant.alleles = {out.flat->reference};
    record.variant.format = {"GT", "PS"};
    record.variant.samples["sample"]["GT"] = {site_gt};
    record.variant.samples["sample"]["PS"] = {"7"};
    record.unflattened_position = 101;
    record.alleles_merged = false;
    for (const FlatStrand& strand : out.flat->strands) {
        record.trav_to_allele[strand.trav] = 0;
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

    SECTION("Differences of the two strands that share reference bases are one record") {
        // Strand 0 leaves out node 4 at the inner site. Strand 1 takes node 6 at the middle
        // site, in place of nodes 3 to 5.
        fixture.phase(fixture.top_key, 0, 0);
        fixture.phase(fixture.middle_key, 0, 1);
        fixture.phase(fixture.inner_key, 1, -1, 1, 0);
        fixture.index();
        const TopSiteFlatRecords top = flatten_top_site(fixture, writer, "0|0");
        REQUIRE(top.flat);
        REQUIRE(top.written == 1);
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat, top.records, 101, "sample"));
        const vcflib::Variant& line = top.records[0];
        REQUIRE(line.position == 105);
        REQUIRE(line.ref == "GCTGAGCA");
        REQUIRE(line.alt == vector<string>{"GGCA", "T"});
        REQUIRE(line.samples.at("sample").at("GT") == vector<string>{"1|2"});
        REQUIRE(line.samples.at("sample").at("PS") == vector<string>{"7"});
        REQUIRE(line.id == ">1>8_0");
        REQUIRE(line.info.at("SB") == vector<string>{"0", "1"});
        REQUIRE(line.infoFlags.at("FL"));
        // No walk spells a flat allele.
        REQUIRE(line.info.at("AT") == vector<string>{".", ".", "."});

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
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat, top.records, 101, "sample"));
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
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat, top.records, 101, "sample"));
        REQUIRE(top.records[0].position == 107);
        REQUIRE(top.records[0].ref == "T");
        REQUIRE(top.records[0].alt == vector<string>{"A"});
        REQUIRE(top.records[0].samples.at("sample").at("GT") == vector<string>{"0|1"});

        // The check fails for records that spell something else.
        vector<vcflib::Variant> wrong = top.records;
        wrong[0].samples["sample"]["GT"] = {"1|0"};
        REQUIRE(!BlockRecordWriter::flat_records_spell(*top.flat, wrong, 101, "sample"));
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
        REQUIRE(BlockRecordWriter::flat_records_spell(*top.flat, top.records, 101, "sample"));
        // The deletion covers the changed base, so the two are one record.
        REQUIRE(top.written == 1);
        REQUIRE(top.records[0].samples.at("sample").at("GT") == vector<string>{"1/2"});
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

TEST_CASE("Flat records of two strands with differences apart are phased by strand",
          "[flat_records]") {
    BlockRecordWriter writer;
    // Strand 0 has a mismatch and an insertion in a repeat. Strand 1 has the same mismatch and a
    // deletion further on.
    const string reference = "GATTACAGGCACACACATTGCCGGAATT";
    const string first = "GATTCCAGGCACACACACATTGCCGGAATT";
    const string second = "GATTCCAGGCACACACATTGCAATT";
    optional<FlatSite> flat = writer.align_flat(reference, {{1, first}, {2, second}});
    REQUIRE(flat);
    SiteRecord record{};
    record.variant.sequenceName = "ref";
    record.variant.position = 1;
    record.variant.id = ">1>2";
    record.variant.ref = reference;
    record.variant.alt = {first, second};
    record.variant.alleles = {reference, first, second};
    record.variant.format = {"GT"};
    record.variant.samples["sample"]["GT"] = {"1|2"};
    record.unflattened_position = 1;
    record.alleles_merged = false;
    record.trav_to_allele = {{0, 0}, {1, 1}, {2, 2}};
    bdsg::HashGraph graph;
    bdsg::PositionOverlay positioned(&graph);
    vector<vcflib::Variant> records;
    REQUIRE(writer.write_flat(positioned, *flat, "sample", nullptr, record, GLLayout::IMajor,
                              [&](vcflib::Variant& line, size_t) {
                                  records.push_back(line);
                                  return true;
                              })
            == 3);
    REQUIRE(BlockRecordWriter::flat_records_spell(*flat, records, 1, "sample"));
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

}
}
