/**
 * \file dinucleotide_machine.cpp
 *
 * Implements DinucleotideMachine
 *
 */

#include "dinucleotide_machine.hpp"

//#define debug_machine

#ifdef debug_machine
#include <bitset>
#endif

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
    // handle transitions to the XN state (the lookup value for N indexes
    // past the rest of the transition table)
    for (size_t i = 128; i < 256; ++i) {
        transition_table[i] = init_state();
    }
    
#ifdef debug_machine
    cerr << "constructed transition table:" << endl;
    for (size_t i = 0; i < 32; ++i) {
        cerr << i << "\t" << bitset<8>(i);
        for (size_t j = 0; j < 4; ++j) {
            cerr << "\t" << bitset<32>(transition_table[4 * i + j]);
        }
        cerr << endl;
    }
#endif
    
    // build a translation table for ASCII nucleotides
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
                // this will cause us to index past the entire table into
                // the XN state
                nt_table[i] = 128;
                break;
        }
    }
}

uint32_t DinucleotideMachine::init_state() const {
    // start in the XN state
    return 1 << 20;
}

uint32_t DinucleotideMachine::update_state(uint32_t state, char next) const {
    // merge the dinucleotide set according to it's final base from positions [15,0]
    uint32_t transition_row = state | (state >> 4);
    transition_row |= (transition_row >> 8);
    // merge in the XN and NA...NT states from positions [20,16]
    transition_row = (transition_row & 0xf) | (state >> 16);
    // do the transitions
    return transition_table[(transition_row << 2) | nt_table[next]];
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
