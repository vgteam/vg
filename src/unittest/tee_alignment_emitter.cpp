/// \file tee_alignment_emitter.cpp
///
/// Unit tests for TeeAlignmentEmitter.

#include <memory>
#include <string>
#include <vector>

#include "catch.hpp"
#include "../tee_alignment_emitter.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A log of what an emitter was given: one line per alignment, naming the
/// method that got it, and one line per extra message.
using EmitterLog = vector<string>;

/// An emitter that writes what it is given into a log it shares with the test,
/// and then clears the alignments and data it was given. If two destinations
/// shared one batch, the second would log empty names.
class LoggingEmitter : public vg::io::AlignmentEmitter {
public:
    explicit LoggingEmitter(shared_ptr<EmitterLog> log) : log(std::move(log)) {}

    void emit_singles(vector<Alignment>&& aln_batch) override {
        for (auto& aln : aln_batch) {
            log->push_back("single " + aln.name());
            aln.Clear();
        }
    }
    void emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) override {
        for (auto& alns : alns_batch) {
            for (auto& aln : alns) {
                log->push_back("mapped single " + aln.name());
                aln.Clear();
            }
        }
    }
    void emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch,
                    vector<int64_t>&& tlen_limit_batch) override {
        for (size_t i = 0; i < aln1_batch.size(); i++) {
            log->push_back("pair " + aln1_batch[i].name() + " " + aln2_batch[i].name() +
                           " " + to_string(tlen_limit_batch.at(i)));
            aln1_batch[i].Clear();
            aln2_batch[i].Clear();
        }
        tlen_limit_batch.clear();
    }
    void emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch,
                           vector<int64_t>&& tlen_limit_batch) override {
        for (size_t i = 0; i < alns1_batch.size(); i++) {
            for (size_t j = 0; j < alns1_batch[i].size(); j++) {
                log->push_back("mapped pair " + alns1_batch[i][j].name() + " " + alns2_batch[i][j].name() +
                               " " + to_string(tlen_limit_batch.at(i)));
                alns1_batch[i][j].Clear();
                alns2_batch[i][j].Clear();
            }
        }
        tlen_limit_batch.clear();
    }
    void emit_extra_message(const string& tag, string&& data) override {
        log->push_back("extra " + tag + " " + data);
        data.clear();
    }

private:
    shared_ptr<EmitterLog> log;
};

/// Make an alignment with the given read name.
Alignment named(const string& name) {
    Alignment aln;
    aln.set_name(name);
    return aln;
}

/// Make a TeeAlignmentEmitter with the given number of LoggingEmitters as
/// destinations, and give it one batch of each kind. Return the log of each
/// destination, in order.
vector<EmitterLog> emit_through_tee(size_t destination_count) {
    vector<shared_ptr<EmitterLog>> logs;
    vector<unique_ptr<vg::io::AlignmentEmitter>> destinations;
    for (size_t i = 0; i < destination_count; i++) {
        logs.push_back(make_shared<EmitterLog>());
        destinations.push_back(make_unique<LoggingEmitter>(logs.back()));
    }

    TeeAlignmentEmitter tee(std::move(destinations));
    tee.emit_extra_message("PARAMS_JSON", "{}");
    tee.emit_singles({named("a"), named("b")});
    tee.emit_mapped_singles({{named("c"), named("c2")}});
    tee.emit_pairs({named("d1")}, {named("d2")}, {500});
    tee.emit_mapped_pairs({{named("e1"), named("f1")}}, {{named("e2"), named("f2")}}, {600});

    vector<EmitterLog> results;
    for (auto& log : logs) {
        results.push_back(*log);
    }
    return results;
}

/// The log lines a LoggingEmitter writes for one extra message and one batch of
/// each kind, holding reads named a through f2.
const EmitterLog EXPECTED_LOG {
    "extra PARAMS_JSON {}",
    "single a", "single b",
    "mapped single c", "mapped single c2",
    "pair d1 d2 500",
    "mapped pair e1 e2 600", "mapped pair f1 f2 600"
};

}

TEST_CASE("TeeAlignmentEmitter gives every destination the same alignments", "[alignment][emitter]") {
    SECTION("One destination gets everything") {
        vector<EmitterLog> logs = emit_through_tee(1);
        REQUIRE(logs.size() == 1);
        REQUIRE(logs[0] == EXPECTED_LOG);
    }

    SECTION("Each of three destinations gets everything, though each clears what it gets") {
        vector<EmitterLog> logs = emit_through_tee(3);
        REQUIRE(logs.size() == 3);
        for (auto& log : logs) {
            REQUIRE(log == EXPECTED_LOG);
        }
    }

    SECTION("With no destinations, emitting does nothing") {
        vector<EmitterLog> logs = emit_through_tee(0);
        REQUIRE(logs.empty());
    }
}

TEST_CASE("TeeAlignmentEmitter refuses a null destination", "[alignment][emitter]") {
    auto log = make_shared<EmitterLog>();
    vector<unique_ptr<vg::io::AlignmentEmitter>> destinations;
    destinations.push_back(make_unique<LoggingEmitter>(log));
    destinations.push_back(nullptr);
    REQUIRE_THROWS_AS(TeeAlignmentEmitter(std::move(destinations)), std::invalid_argument);
}

}
}
