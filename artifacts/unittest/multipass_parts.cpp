/// \file multipass_parts.cpp
/// Tests for the parts of the read-likelihood caller that sort buffered VCF records, give records
/// their nesting tags, and phase reads and nested sites.
///
/// Records sharing a position are rare in real output, so the order of ties is checked here.

#include <algorithm>
#include <cmath>
#include <sstream>
#include <map>
#include <functional>
#include <string>
#include <string_view>
#include <vector>


#include "../../src/unittest/catch.hpp"
#include "support/decomposition_pair.hpp"
#include "../edit_snarl_decomposition_overlay.hpp"
#include "../../src/integrated_snarl_finder.hpp"
#include "support/nested_graph.hpp"
#include "../../src/graph_caller.hpp"
#include "../../src/flow_caller.hpp"
#include "../../src/sorted_vcf_writer.hpp"
#include "../descent_policy.hpp"
#include "../phase_table.hpp"
#include "../read_strand_table.hpp"
#include "../staged_site.hpp"
#include "../vcf_record.hpp"
#include "../../src/utility.hpp"

namespace vg {
namespace unittest {

using namespace vg::multipass;

using namespace std;

TEST_CASE("The buffered record order is total, so two runs cannot disagree", "[graph_caller]") {
    using BufferedRecordKey = SortedVCFWriter::BufferedRecordKey;

    // Two records at the same position, distinguished only by the snarl they came from. This is the
    // case (contig, POS) cannot order, and it is not rare: a nested site sits at or near its
    // parent's position, and flattening can put two snarls on one anchor base.
    const BufferedRecordKey a{"chr20", 1000, ">1>5"};
    const BufferedRecordKey b{"chr20", 1000, ">6>9"};

    // Antisymmetry: without the ID in the key, `<` would be false in both directions for this
    // pair, so `std::sort` could leave them in input order.
    REQUIRE(a < b);
    REQUIRE_FALSE(b < a);

    // Irreflexivity, since a strict weak ordering needs it and a `<=` typo here would make
    // std::sort's behaviour undefined rather than merely wrong.
    REQUIRE_FALSE(a < a);

    SECTION("the earlier field always dominates the later one") {
        // Position beats ID: a later ID at an earlier position still sorts first, so the tie-break
        // cannot reorder the file.
        REQUIRE(BufferedRecordKey{"chr20", 999, ">9>9"} < BufferedRecordKey{"chr20", 1000, ">1>1"});
        // Contig beats both.
        REQUIRE(BufferedRecordKey{"chr1", 5000, ">9>9"} < BufferedRecordKey{"chr20", 1, ">1>1"});
    }

    SECTION("every input permutation sorts to one output") {
        // Sorting each permutation of a set containing a tie must give the same sequence every
        // time.
        vector<BufferedRecordKey> keys{
            {"chr20", 1000, ">6>9"},
            {"chr20", 1000, ">1>5"},
            {"chr20", 900, ">2>3"},
            {"chr21", 10, ">4>7"},
        };
        sort(keys.begin(), keys.end());
        const vector<string> want{">2>3", ">1>5", ">6>9", ">4>7"};

        // Start from the first permutation in id order, so next_permutation walks all of them
        // rather than reporting exhaustion on its first call.
        vector<BufferedRecordKey> permuted = keys;
        sort(permuted.begin(), permuted.end(),
             [](const BufferedRecordKey& x, const BufferedRecordKey& y) { return x.id < y.id; });
        size_t checked = 0;
        do {
            vector<BufferedRecordKey> copy = permuted;
            sort(copy.begin(), copy.end());
            for (size_t i = 0; i < want.size(); ++i) {
                REQUIRE(copy[i].id == want[i]);
            }
            ++checked;
        } while (next_permutation(permuted.begin(), permuted.end(),
                                  [](const BufferedRecordKey& x, const BufferedRecordKey& y) {
                                      return x.id < y.id;
                                  }));
        // All 24 permutations of four distinct records, so the claim is exhaustive rather than
        // sampled.
        REQUIRE(checked == 24);
    }
}

namespace {

/// Return a record on contig "chr" at `position` with ID `id` and REF `ref`, called 0/1.
vcflib::Variant nesting_record(const string& id, int64_t position, const string& ref) {
    vcflib::Variant record;
    record.sequenceName = "chr";
    record.position = position;
    record.id = id;
    record.ref = ref;
    record.alt = {"A"};
    record.alleles = {ref, "A"};
    record.filter = "PASS";
    record.info["AT"] = {">1", ">2"};
    record.format = {"GT"};
    record.samples["sample"]["GT"] = {"0/1"};
    return record;
}

/// Return the INFO fields of each record line of a VCF, by the record's ID.
map<string, map<string, string>> info_by_id(const string& vcf) {
    map<string, map<string, string>> fields;
    stringstream lines(vcf);
    string line;
    while (getline(lines, line)) {
        if (line.starts_with("#")) {
            continue;
        }
        vector<string> columns = split_delims(line, "\t");
        for (const string& field : split_delims(columns.at(7), ";")) {
            size_t equals = field.find('=');
            fields[columns.at(2)][field.substr(0, equals)] =
                equals == string::npos ? "" : field.substr(equals + 1);
        }
    }
    return fields;
}

}

TEST_CASE("The nesting tags place each record under the records of the sites around it",
          "[graph_caller]") {
    NestedGraph nested;
    const HandleGraph& graph = nested.graph;
    IntegratedSnarlFinder finder(nested.graph);
    DecompositionPair pair(nested.graph, finder);
    EditSnarlDecompositionOverlay overlay(pair.adapter, nested.graph);
    // An edit site inside node 4, which lies inside 3..5, inside 2..7, inside 1..8.
    const net_handle_t edit = overlay.add_edit(nested.handle(4), 1, nested.handle(4), 2);
    auto record_id_of_site = [&](const net_handle_t& site) {
        if (overlay.is_edit(site)) {
            return string("4:1:G");
        }
        return site_name(graph, overlay.get_site_start(site), overlay.get_site_end(site), nullptr,
                         false);
    };
    REQUIRE(record_id_of_site(edit) == "4:1:G");

    SortedVCFWriter output("sample");
    output.set_nested(true);
    // The top-level site, the middle site written read backward, and the edit site. The site
    // 3..5 between the edit and the middle site has no record.
    vcflib::Variant top = nesting_record(">1>8", 3, "ATGCTGAGCATG");
    vcflib::Variant middle = nesting_record("<7<2", 4, "TGCTGAG");
    vcflib::Variant inner = nesting_record("4:1:G", 7, "T");
    REQUIRE(output.add_variant(top));
    REQUIRE(output.add_variant(middle));
    REQUIRE(output.add_variant(inner));

    stringstream vcf;
    output.write_variants(vcf, &overlay, &graph, record_id_of_site);
    map<string, map<string, string>> info = info_by_id(vcf.str());
    REQUIRE(info.size() == 3);

    // The top-level record has no parent, and is its own place on the reference.
    REQUIRE(info[">1>8"]["LV"] == "0");
    REQUIRE(info[">1>8"]["CH"] == "0");
    REQUIRE_FALSE(info[">1>8"].contains("PS"));
    REQUIRE(info[">1>8"]["RC"] == "chr");
    REQUIRE(info[">1>8"]["RS"] == "3");
    REQUIRE(info[">1>8"]["RD"] == "15");

    // The middle record is found though it is written read backward.
    REQUIRE(info["<7<2"]["LV"] == "1");
    REQUIRE(info["<7<2"]["CH"] == "0");
    REQUIRE(info["<7<2"]["PS"] == ">1>8");
    REQUIRE(info["<7<2"]["RC"] == "chr");
    REQUIRE(info["<7<2"]["RS"] == "3");
    REQUIRE(info["<7<2"]["RD"] == "15");

    // The edit's record is found by its own ID. Its parent is the nearest site around it with a
    // record, and it is placed on the reference where the top-level record is.
    REQUIRE(info["4:1:G"]["LV"] == "2");
    REQUIRE(info["4:1:G"]["CH"] == "0");
    REQUIRE(info["4:1:G"]["PS"] == "<7<2");
    REQUIRE(info["4:1:G"]["RC"] == "chr");
    REQUIRE(info["4:1:G"]["RS"] == "3");
    REQUIRE(info["4:1:G"]["RD"] == "15");
}

TEST_CASE("A strand from another phase set is NaN in the anchor path, and no strand is 0",
          "[graph_caller]") {
    // When homozygous sites are split for the anchors, a NaN read is dropped from a split site and
    // a 0 read is placed by a coin, so the two must stay distinct here. Re-genotyping does not use
    // this lookup, and gives both 0.
    auto key = [](const string& name) { return (uint64_t)std::hash<std::string_view>{}(name); };
    // Each read carries the allele in slot 0, so each site it is on points it to strand 0.
    auto site = [&](size_t record_key, size_t phase_set, const vector<string>& reads) {
        PhaseSite s;
        s.record_key = record_key;
        s.phase_set = phase_set;
        for (const string& name : reads) {
            s.read_key.push_back(key(name));
            s.q0.push_back(0.9f);
            s.c.push_back(1.0f);
        }
        return s;
    };
    ReadStrandTable strands;
    // "here" is on two sites of phase set 7, "elsewhere" on one site of phase set 9, and "both" on
    // one site of each.
    strands.sites() = {site(10, 7, {"here", "both"}), site(11, 7, {"here"}),
                       site(12, 9, {"elsewhere", "both"})};
    // Site 1 is in phase set 7 and has no read evidence of its own. Site 2 has no phase set.
    auto call = [](size_t record_key, size_t phase_set) {
        ReadStrandTable::PhaseCall pc;
        pc.record_key = record_key;
        pc.phase_set = phase_set;
        return pc;
    };
    strands.build_lambda({call(1, 7), call(12, 9)},
                         [](const string&, size_t phase_set) { return phase_set; },
                         StrandCalibration{1.0, 1.0}, RegenotypeParams());

    const double here = strands.read_strand_log_odds(1, "here");
    REQUIRE(std::isfinite(here));
    REQUIRE(here > 0.0);
    REQUIRE(std::isnan(strands.read_strand_log_odds(1, "elsewhere")));
    REQUIRE(std::isnan(strands.read_strand_log_odds(1, "both")));
    REQUIRE(strands.read_strand_log_odds(1, "unseen") == 0.0);

    // A site with no phase set takes any read found in one phase set, but not one found in two.
    REQUIRE(strands.read_strand_log_odds(2, "elsewhere") > 0.0);
    REQUIRE(std::isnan(strands.read_strand_log_odds(2, "both")));

    // A site's own evidence is left out, so a read on no other site has no strand there.
    REQUIRE(strands.read_strand_log_odds(10, "here") > 0.0);
    REQUIRE(strands.read_strand_log_odds(10, "here") < here);
    REQUIRE(strands.read_strand_log_odds(12, "elsewhere") == 0.0);
}

TEST_CASE("A parent's phase swap carries its nested strands and their haplotypes with it",
          "[graph_caller]") {
    const size_t W = LinkageModel::NO_HAPLOTYPE;
    // Diploid parent 1, swapped by read phasing. Its ploidy-1 child 2 is on strand 1 with
    // haplotype 7, and 2's own ploidy-1 child 3 is on strand 0 with haplotype 5. Diploid child 4 is
    // not swapped, so its ploidy-1 child 5 keeps its strand.
    auto call = [](size_t key, size_t ploidy, int strand, size_t first, size_t second) {
        PhaseTable::PhaseCall pc;
        pc.record_key = key;
        pc.ploidy = ploidy;
        pc.nested_strand = (int8_t)strand;
        pc.hap_first = first;
        pc.hap_second = second;
        return pc;
    };
    PhaseTable table;
    vector<PhaseTable::PhaseCall>& phased = table.calls();
    phased = {
        call(1, 2, -1, 3, 4), call(2, 1, 1, W, 7), call(3, 1, 0, 5, W),
        call(4, 2, -1, 3, 4), call(5, 1, 0, 6, W)};
    // Out of level order, as the staged sites can be.
    vector<PhaseTable::NestedLink> links = {
        {3, 2, 2}, {5, 4, 2}, {2, 1, 1}, {4, 1, 1}, {1, 0, 0}};
    const unordered_set<size_t> flips = {1};

    SECTION("each strand moves, and its haplotype moves to the slot the strand names") {
        REQUIRE(table.swap_strands(flips, links) == 2);
        // The swapped parent's own haplotypes change strands.
        REQUIRE(phased[0].hap_first == 4);
        REQUIRE(phased[0].hap_second == 3);
        REQUIRE(phased[1].nested_strand == 0);
        REQUIRE(phased[1].hap_first == 7);
        REQUIRE(phased[1].hap_second == W);
        REQUIRE(phased[2].nested_strand == 1);
        REQUIRE(phased[2].hap_first == W);
        REQUIRE(phased[2].hap_second == 5);
        REQUIRE(phased[4].nested_strand == 0);
        REQUIRE(phased[4].hap_first == 6);
    }
    SECTION("a chain left out of the links stops the swap reaching its children") {
        links.erase(links.begin() + 2);   // chain 2
        REQUIRE(table.swap_strands(flips, links) == 0);
        REQUIRE(phased[2].nested_strand == 0);
    }
}

}
}
