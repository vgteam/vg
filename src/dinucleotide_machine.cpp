/**
 * \file dinucleotide_machine.cpp
 *
 * Implements DinucleotideMachine
 *
 */

#include "dinucleotide_machine.hpp"

namespace vg {
    
DinucleotideMachine::DinucleotideMachine() {
    // build the transition table for the state machine
    for (size_t i = 0; i < 16; ++i) {
        uint16_t base = 0;
        for (size_t j = 0; j < 4; ++j) {
            if (i & (1 << j)) {
                base |= 1 << (4 * j);
            }
        }
        for (size_t j = 0; j < 4; ++j) {
            // the dinucleotide set ending in j
            transition_table[4 * i + j] = base << j;
            // the dinucleotide set ending in j including Nj
            transition_table[4 * i + j + 64] = (base << j) | (1 << (16 + j));
        }
    }
    
    // build a translation table for
    for (size_t i = 0; i < 256; ++i) {
        switch (i) {
            case 'a':
            case 'A':
                nt_table[i] = 0;
                break;
            case 'c':
            case 'C':
                nt_table[i] = 1;
                break;
            case 'g':
            case 'G':
                nt_table[i] = 2;
                break;
            case 't':
            case 'T':
                nt_table[i] = 3;
                break;
            default:
                nt_table[i] = numeric_limits<uint32_t>::max();
                break;
        }
    }
}

uint32_t DinucleotideMachine::init_state() const {
    // start in the XN state
    return 1 << 20;
}

uint32_t DinucleotideMachine::update_state(uint32_t state, char next) const {
    if (next == 'N') {
        // return the XN state
        return 1 << 20;
    }
    else {
        // merge the dinucleotide set according to it's final base
        uint32_t next_state = state | (state >> 1);
        next_state |= (next_state >> 2);
        // merge in the XN and NA...NT states and transition
        return transition_table[(((next_state & 0xff) | (state >> 16)) << 2) | nt_table[next]];
    }
}


uint32_t DinucleotideMachine::merge_state(uint32_t state_1, uint32_t state_2) const {
    return state_1 | state_2;
}

bool DinucleotideMachine::matches(uint32_t state, const char* dinucleotide) const {
    return state & (1 << ((nt_table[dinucleotide[0]] << 2) | nt_table[dinucleotide[1]]));
}

bool DinucleotideMachine::matches(uint32_t state, const string& dinucleotide) const {
    return matches(state, dinucleotide.c_str());
}

}
