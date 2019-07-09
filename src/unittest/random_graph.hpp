#include "handle.hpp"
#include <random>
#include <time.h>

namespace vg{
namespace unittest{

/// Create a random graph by adding variation to a sequence of length seq_size
/// variant_len is the mean length of a larger variation and variant_count
/// is the number of variations added to the graph
void random_graph(int64_t seq_size, int64_t variant_len, int64_t variant_count,
                  MutablePathMutableHandleGraph* graph);

}
}
