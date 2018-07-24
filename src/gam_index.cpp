#include "gam_index.hpp"

#include <iostream>

namespace vg {

using namespace std;

auto GAMIndex::bins_of_id(id_t id) -> vector<bin_t> {
    vector<bin_t> to_return;
    
    // How many bits are in the number? We get a bin per bit.
    auto number_bits = CHAR_BIT * sizeof(bin_t);
    
    // We need to keep in mind that shifting *all* the bits out of a number in
    // one go is undefined behavior.
    
    // Bin number consists of an index which is a prefix of the number
    bin_t bin_index = ((bin_t)id) >> 1;
    // And an offset to keep from colliding with other bins.
    bin_t bin_offset = ((bin_t)~0) >> 1;
    
    for (int i = 0; i < number_bits; i++) {
        
#ifdef debug
        cerr << hex << id << " " << i << " " << bin_index << " " << bin_offset << " " << (bin_index + bin_offset) << endl;
#endif
    
        to_return.push_back(bin_index + bin_offset);
        bin_index = bin_index >> 1;
        bin_offset = bin_offset >> 1;
    }
    
    return to_return;
}
    
auto GAMIndex::common_bin(id_t a, id_t b) -> bin_t {
    // Convert to unsigned numbers
    bin_t a_bin = a;
    bin_t b_bin = b;
    
    // Define the offset for the bin
    bin_t offset = ((bin_t)~0);
    
    // We're just going to pop off bits until we find the common prefix.
    // Always pop off one bit, even if we are binning a number and itself.
    // TODO: Find a faster way to do this with the appropriate instruction intrinsics.
    do {
        a_bin = a_bin >> 1;
        b_bin = b_bin >> 1;
        offset = offset >> 1;
    } while(a_bin != b_bin);
    return a_bin + offset;
}

auto GAMIndex::window_of_id(id_t id) -> window_t  {
    return id >> WINDOW_SHIFT;
}


}
