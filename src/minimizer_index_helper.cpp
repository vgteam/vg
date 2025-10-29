/**
 * \file minimizer_index_helper.cpp
 * Implementation of shared logic for minimizer index construction.
 */

#include "minimizer_index_helper.hpp"

namespace mi_helper {

size_t trailing_zeros(size_t value) {
  size_t result = 0;
  if (value == 0) {
    return result;
  }
  while ((value & 1) == 0) {
    value >>= 1;
    result++;
  }
  return result;
}

size_t estimate_hash_table_size(const gbwtgraph::GBZ &gbz, bool progress) {
  if (progress) {
    std::cerr << "Estimating genome size" << std::endl;
  }
  size_t genome_size = 0;

  if (gbz.graph.get_path_count() > 0) {
    gbz.graph.for_each_path_handle([&](const path_handle_t &path_handle) {
      gbz.graph.for_each_step_in_path(
          path_handle, [&](const step_handle_t &step_handle) {
            handle_t handle = gbz.graph.get_handle_of_step(step_handle);
            genome_size += gbz.graph.get_length(handle);
          });
    });
    if (progress) {
      std::cerr << "Estimated size based on reference / generic paths: "
                << genome_size << std::endl;
    }
  }

  if (genome_size == 0) {
    gbz.graph.for_each_handle([&](const handle_t &handle) {
      genome_size += gbz.graph.get_length(handle);
    });
    if (progress) {
      std::cerr << "Estimated size based on total sequence length: "
                << genome_size << std::endl;
    }
  }

  // Genome size / 2 should be a reasonably tight upper bound for the number of
  // kmers with any specific base in the middle position.
  size_t hash_table_size =
      gbwtgraph::KmerIndex<gbwtgraph::Key64, gbwtgraph::Position>::minimum_size(
          genome_size / 2);
  if (progress) {
    std::cerr << "Estimated hash table size: 2^"
              << trailing_zeros(hash_table_size) << std::endl;
  }

  return hash_table_size;
}

} // namespace mi_helper