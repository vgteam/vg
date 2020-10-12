/**
 * \file dinucleotide_machine.hpp
 *
 * Defines a nondeterministic finite automaton over dinucleotides
 *
 */
#ifndef VG_DINUCLEOTIDE_MACHINE_GRAPH_HPP_INCLUDED
#define VG_DINUCLEOTIDE_MACHINE_GRAPH_HPP_INCLUDED

#include <string>
#include <limits>
#include <cstdint>
#include <iostream>

namespace vg {

using namespace std;

/*
 * Represents a non-deterministic finite automaton whose states
 * correspond to dinucleotides
 */
class DinucleotideMachine {
public:
    DinucleotideMachine();
    ~DinucleotideMachine() = default;
    
    /// Return an empty dinucleotide set
    uint32_t init_state() const;
    
    /// Get the dinucleotide set that results from extending the set by the given character
    /// Ns are valid, but will result in no matches
    uint32_t update_state(uint32_t state, char next) const;
    
    /// Get the union of two dinucleotide sets
    uint32_t merge_state(uint32_t state_1, uint32_t state_2) const;
    
    /// Return true if the set includes this dinucleotide. Only valid for dinucleotides of ACGT (never N).
    bool matches(uint32_t state, const char* dinucleotide) const;
    
    /// Same semantics as above
    bool matches(uint32_t state, const string& dinucleotide) const;
    
private:
    
    // lookup table for transitions
    uint32_t transition_table[256];
    // ASCII-indexed conversion from char to table index
    uint32_t nt_table[256];
};

}

#endif
