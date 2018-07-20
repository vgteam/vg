#include "gam_index.hpp"

namespace vg {

using namespace std;

auto GAMIndex::bins_of_id(id_t id) -> vector<bin_t> {
    return vector<bin_t>();
}
    
auto GAMIndex::common_bin(id_t a, id_t b) -> bin_t {
    return 0;
}

auto GAMIndex::window_of_id(id_t id) -> window_t  {
    return 0;
}


}
