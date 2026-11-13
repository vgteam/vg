/// \file tee_alignment_emitter.cpp
///
/// Unit tests for TeeAlignmentEmitter.

#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "catch.hpp"
#include "../tee_alignment_emitter.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// Test AlignmentEmitter that makes a text log of the alginments it got.
/// Makes sure to clear anything it can move out of, to annoy the source.
/// Must outlive its logging destination.
class LoggingEmitter : public vg::io::AlignmentEmitter {
public:
    explicit LoggingEmitter(std::ostream& log) : log(log) {
        // Nothing to do!
    }

    void emit_singles(vector<Alignment>&& aln_batch) override {
        for (auto& aln : aln_batch) {
            log << "single " << aln.name() << std::endl;
            aln.Clear();
        }
    }

    void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) override {
        for (auto& alns : alns_batch) {
            for (auto& aln : alns) {
                log << "mapped single " << aln.name() << std::endl;
                aln.Clear();
            }
        }
    }

    void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
                    vector<int64_t>&& tlen_limit_batch) override {
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            log << "pair " << aln1_batch[i].name() << " " 
                << aln2_batch[i].name() << " "
                << tlen_limit_batch.at(i) << std::endl;
            aln1_batch[i].Clear();
            aln2_batch[i].Clear();
        }
        tlen_limit_batch.clear();
    }

    void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch,
                           vector<int64_t>&& tlen_limit_batch) override {
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            for (size_t j = 0; j < alns1_batch[i].size(); j++) {
                log << "mapped pair " << alns1_batch[i][j].name() << " "
                    << alns2_batch[i][j].name() << " "
                    << tlen_limit_batch.at(i) << std::endl;
                alns1_batch[i][j].Clear();
                alns2_batch[i][j].Clear();
            }
        }
        tlen_limit_batch.clear();
    }

    void emit_extra_message(const string& tag, string&& data) override {
        log << "extra " << tag << " " << data << std::endl;
        data.clear();
    }

private:
    std::ostream& log;
};

/// Make an alignment with the given read name.
Alignment named(const string& name) {
    Alignment aln;
    aln.set_name(name);
    return aln;
}

/// Make a TeeAlignmentEmitter with the given number of LoggingEmitters as
/// destinations, and exercise all its functions.
///
/// Return all the logs, which should match.
vector<std::string> emit_through_tee(size_t destination_count) {
    vector<std::stringstream> logs;
    vector<unique_ptr<vg::io::AlignmentEmitter>> destinations;
    logs.resize(destination_count);
    for (size_t i = 0; i < destination_count; i++) {
        destinations.push_back(make_unique<LoggingEmitter>(logs[i]));
    }

    TeeAlignmentEmitter tee(std::move(destinations));
    tee.emit_extra_message("PARAMS_JSON", "{}");
    tee.emit_singles({named("a"), named("b")});
    tee.emit_mapped_singles({{named("c"), named("c2")}});
    tee.emit_pairs({named("d1")}, {named("d2")}, {500});
    tee.emit_mapped_pairs({{named("e1"), named("f1")}}, {{named("e2"), named("f2")}}, {600});

    vector<std::string> results;
    for (auto& log : logs) {
        results.push_back(log.str());

    }
    return results;
}

/// The expected log for each emitter from emit_through_tee()
const std::string EXPECTED_LOG = R"(extra PARAMS_JSON {}
single a
single b
mapped single c
mapped single c2
pair d1 d2 500
mapped pair e1 e2 600
mapped pair f1 f2 600
)";

}

TEST_CASE("TeeAlignmentEmitter gives every destination the same alignments", "[alignment][emitter]") {
    SECTION("Each of three destinations gets everything, though each clears what it gets") {
        vector<std::string> logs = emit_through_tee(3);
        REQUIRE(logs.size() == 3);
        for (auto& log : logs) {
            REQUIRE(log == EXPECTED_LOG);
        }
    }

    SECTION("With no destinations, emitting is allowed") {
        vector<std::string> logs = emit_through_tee(0);
        REQUIRE(logs.empty());
    }
}

TEST_CASE("TeeAlignmentEmitter refuses a null destination", "[alignment][emitter]") {
    vector<unique_ptr<vg::io::AlignmentEmitter>> destinations;
    destinations.push_back(nullptr);
    REQUIRE_THROWS_AS(TeeAlignmentEmitter(std::move(destinations)), std::invalid_argument);
}

}
}
