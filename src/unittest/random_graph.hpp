#include "vg.hpp"
#include <random>
#include <time.h>

namespace vg{
namespace unittest{

/// Create a random graph for a sequence of length seqSize
/// variantLen is the mean length of a larger variation and variationCount
/// is the number of variations in the graph
VG randomGraph(int64_t seqSize, int64_t variantLen, int64_t variantCount);

}
}
