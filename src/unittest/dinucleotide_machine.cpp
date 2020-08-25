/// \file dinucleotide_machine.cpp
///  
/// Unit tests for the DinucleotideMachine
///

#include <iostream>
#include <random>

#include "../dinucleotide_machine.hpp"
#include "../statistics.hpp"
#include "catch.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {
using namespace std;

using bdsg::HashGraph;

TEST_CASE("DinucleotideMachine can correctly detect simple dinucleotides",
          "[dinucleotide]") {
    
    DinucleotideMachine machine;
    
    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            auto state = machine.init_state();
            
            state = machine.update_state(state, nt1);
            state = machine.update_state(state, nt2);
            
            for (char test_nt1 : string("ACGT")) {
                for (char test_nt2 : string("ACGT")) {
                    string di_nt{test_nt1, test_nt2};
                    if (test_nt1 == nt1 && test_nt2 == nt2) {
                        REQUIRE(machine.matches(state, di_nt));
                    }
                    else{
                        REQUIRE(!machine.matches(state, di_nt));
                    }
                }
            }
        }
    }
}

TEST_CASE("DinucleotideMachine works correctly on random sequences of nucleotides",
          "[dinucleotide]") {
    
    DinucleotideMachine machine;
    
    string nts = "ACGT";
    default_random_engine gen(981727497);
    uniform_int_distribution<int> distr(0, 3);
    
    auto state = machine.init_state();
    
    string di_nt = "AA";
    state = machine.update_state(state, di_nt[0]);
    state = machine.update_state(state, di_nt[1]);
    
    size_t n_trials = 1000;
    for (size_t i = 0; i < n_trials; ++i) {
        di_nt[0] = di_nt[1];
        di_nt[1] = nts[distr(gen)];
        state = machine.update_state(state, di_nt[1]);
        REQUIRE(machine.matches(state, di_nt));
    }
}

TEST_CASE("DinucleotideMachine can merge states correctly", "[dinucleotide]") {
    
    DinucleotideMachine machine;
    
    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            for (char nt3 : string("ACGT")) {
                
                auto state1 = machine.init_state();
                auto state2 = machine.init_state();
                auto state3 = machine.init_state();
                auto state4 = machine.init_state();
                
                state1 = machine.update_state(state1, nt1);
                state2 = machine.update_state(state2, nt2);
                state1 = machine.merge_state(state1, state2);
                state1 = machine.update_state(state1, nt3);
                
                state3 = machine.update_state(state3, nt1);
                state4 = machine.update_state(state3, nt2);
                state3 = machine.update_state(state3, nt3);
                state3 = machine.merge_state(state3, state4);
                
                for (char test_nt1 : string("ACGT")) {
                    for (char test_nt2 : string("ACGT")) {
                        string test_di_nt{test_nt1, test_nt2};
                        if ((test_nt1 == nt1 || test_nt1 == nt2) && test_nt2 == nt3) {
                            REQUIRE(machine.matches(state1, test_di_nt));
                        }
                        else {
                            REQUIRE(!machine.matches(state1, test_di_nt));
                        }
                        if (test_nt1 == nt1 && (test_nt2 == nt2 || test_nt2 == nt3)) {
                            REQUIRE(machine.matches(state3, test_di_nt));
                        }
                        else {
                            REQUIRE(!machine.matches(state3, test_di_nt));
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("DinucleotideMachine handles N correctly", "[dinucleotide]") {
    
    DinucleotideMachine machine;
    
    auto state = machine.init_state();
    state = machine.update_state(state, 'A');
    state = machine.update_state(state, 'A');
    state = machine.update_state(state, 'N');
    
    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            string di_nt{nt1, nt2};
            REQUIRE(!machine.matches(state, di_nt));
        }
    }
    
    state = machine.update_state(state, 'N');

    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            string di_nt{nt1, nt2};
            REQUIRE(!machine.matches(state, di_nt));
        }
    }
    
    state = machine.update_state(state, 'A');
    
    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            string di_nt{nt1, nt2};
            REQUIRE(!machine.matches(state, di_nt));
        }
    }
    
    state = machine.update_state(state, 'A');
    
    for (char nt1 : string("ACGT")) {
        for (char nt2 : string("ACGT")) {
            string di_nt{nt1, nt2};
            if (nt1 == 'A' && nt2 == 'A') {
                REQUIRE(machine.matches(state, di_nt));
            }
            else {
                REQUIRE(!machine.matches(state, di_nt));
            }
        }
    }
    
}

}
}
        
