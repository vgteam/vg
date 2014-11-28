#include "mapper.h"

namespace vg {

Alignment& Mapper::map(Alignment& alignment) {
    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }
    VariantGraph& graph;
    
    return alignment;
}

}
