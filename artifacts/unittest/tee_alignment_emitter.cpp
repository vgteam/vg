/// \file tee_alignment_emitter.cpp
///
/// Unit tests for TeeAlignmentEmitter.

#include <memory>
#include <string>
#include <vector>

#include "../../src/unittest/catch.hpp"
#include "../tee_alignment_emitter.hpp"

namespace vg {
namespace unittest {

using namespace std;

namespace {

/// A log of what an emitter was given: one line per alignment, naming the
/// method that received it, and one line per extra message.
using EmitterLog = vector<string>;

/// An emitter that writes what it is given into a log it shares with the test,
/// and then clears the alignments it was given, as an emitter that moves from
/// its batches would leave them.
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

}

TEST_CASE("TeeAlignmentEmitter gives both emitters the same alignments", "[alignment][emitter]") {
    auto primary_log = make_shared<EmitterLog>();
    auto copy_log = make_shared<EmitterLog>();
    TeeAlignmentEmitter tee(make_unique<LoggingEmitter>(primary_log), make_unique<LoggingEmitter>(copy_log));

    tee.emit_extra_message("PARAMS_JSON", "{}");
    tee.emit_singles({named("a"), named("b")});
    tee.emit_mapped_singles({{named("c"), named("c2")}});
    tee.emit_pairs({named("d1")}, {named("d2")}, {500});
    tee.emit_mapped_pairs({{named("e1"), named("f1")}}, {{named("e2"), named("f2")}}, {600});

    EmitterLog expected {
        "extra PARAMS_JSON {}",
        "single a", "single b",
        "mapped single c", "mapped single c2",
        "pair d1 d2 500",
        "mapped pair e1 e2 600", "mapped pair f1 f2 600"
    };
    REQUIRE(*primary_log == expected);
    REQUIRE(*copy_log == expected);
}

TEST_CASE("TeeAlignmentEmitter needs two emitters", "[alignment][emitter]") {
    auto log = make_shared<EmitterLog>();
    REQUIRE_THROWS_AS(TeeAlignmentEmitter(make_unique<LoggingEmitter>(log), nullptr), std::invalid_argument);
    REQUIRE_THROWS_AS(TeeAlignmentEmitter(nullptr, make_unique<LoggingEmitter>(log)), std::invalid_argument);
}

}
}
