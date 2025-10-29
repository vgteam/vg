/**
\file minimizer_index_helper.hpp
Shared logic for minimizer index construction.
*/
#ifndef VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED
#define VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED

#include <fstream>
#include <functional>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/index.h>
#include <mutex>
#include <omp.h>
#include <string>
#include <type_traits>
#include <vector>

#include "progressive.hpp"
#include "snarl_distance_index.hpp"
#include "vg/io/vpkg.hpp"
#include "zip_code.hpp"

namespace vg {

/**
 * Helper class for building minimizer indexes with zipcodes.
 *
 * Typical usage:
 *   MinimizerIndexHelper helper(gbz, index);
 *   helper.set_frequent_kmers(k, threshold);
 *   helper.build(distance_index, zipcode_filename, output_filename);
 *
 * @tparam IndexType The type of minimizer index to build
 * @tparam PayloadType The type of payload to store in the index
 */
template <typename IndexType, typename PayloadType>
class MinimizerIndexHelper : public Progressive {
public:
  /**
   * Construct a helper for building a minimizer index.
   *
   * @param gbz The GBZ graph to index
   * @param index The minimizer index to populate
   * @param show_progress If true, display progress information
   */
  MinimizerIndexHelper(const gbwtgraph::GBZ *gbz, IndexType &index,
                       bool show_progress = false);

  /**
   * Identify and mark frequent kmers in the index.
   *
   * Finds kmers that occur more than the specified threshold and set them as
   * frequent in the index.
   */
  void set_frequent_kmers(int k, size_t threshold,
                          bool space_efficient_counting = false,
                          size_t hash_table_size = 0, size_t iterations = 1);

  /**
   * Build a minimizer index with type IndexType, holding payloads of type
   * PayloadType.
   *
   * If distance_index is set, uses that distance index to make payloads.
   *
   * If zipcode_name is nonempty, sotres oversize zipcodes to that file.
   *
   * Writes the minimizer index to the file named output_name; the filename must
   * be nonempty.
   *
   * Logs progress if progress is true.
   */
  void build(const bdsg::SnarlDistanceIndex *distance_index,
             const std::string &zipcode_name, const std::string &output_name);

  /**
   * Estimate an appropriate hash table size based on the total length of the
   * sequences in the GBZ graph.
   */
  size_t estimate_hash_table_size() const;

private:
  const gbwtgraph::GBZ *gbz;
  IndexType &index;

  /**
   * Count the number of trailing zero bits in a value.
   */
  static size_t trailing_zeros(size_t value);
};

// Template implementation

template <typename IndexType, typename PayloadType>
MinimizerIndexHelper<IndexType, PayloadType>::MinimizerIndexHelper(
    const gbwtgraph::GBZ *gbz, IndexType &index, bool show_progress)
    : Progressive(), gbz(gbz), index(index) {
  this->show_progress = show_progress;
}

template <typename IndexType, typename PayloadType>
void MinimizerIndexHelper<IndexType, PayloadType>::set_frequent_kmers(
    int k, size_t threshold, bool space_efficient_counting,
    size_t hash_table_size, size_t iterations) {

  std::vector<gbwtgraph::Key64> frequent_kmers;

  if (show_progress) {
    std::string algorithm =
        (space_efficient_counting ? "space-efficient" : "fast");
    std::cerr << "Finding frequent kmers using the " << algorithm
              << " algorithm" << std::endl;
  }

  if (hash_table_size == 0) {
    hash_table_size = estimate_hash_table_size();
  }

  frequent_kmers = gbwtgraph::frequent_kmers<gbwtgraph::Key64>(
      gbz->graph, k, threshold, space_efficient_counting, hash_table_size);

  if (show_progress) {
    std::cerr << "Found " << frequent_kmers.size() << " kmers with more than "
              << threshold << " hits" << std::endl;
  }

  if (!frequent_kmers.empty()) {
    index.add_frequent_kmers(frequent_kmers, iterations);
  }
}

template <typename IndexType, typename PayloadType>
void MinimizerIndexHelper<IndexType, PayloadType>::build(
    const bdsg::SnarlDistanceIndex *distance_index,
    const std::string &zipcode_name, const std::string &output_name) {

  // Zipcodes
  ZipCodeCollection oversized_zipcodes;

  // Map node id to what gets stored in the payload
  hash_map<id_t, PayloadType> node_id_to_payload;
  node_id_to_payload.reserve(gbz->graph.max_node_id() -
                             gbz->graph.min_node_id());

  // Mutexes for protecting shared data structures in parallel regions
  std::mutex payload_mutex;
  std::mutex zipcode_mutex;

  // Build the index
  if (show_progress) {
    std::cerr << "Building MinimizerIndex with k = " << index.k();
    if (index.uses_syncmers()) {
      std::cerr << ", s = " << index.s();
    } else {
      std::cerr << ", w = " << index.w();
    }
    if (index.uses_weighted_minimizers()) {
      std::cerr << ", W";
    }
    std::cerr << " using : ";
    if (std::is_same<PayloadType, gbwtgraph::Payload>::value) {
      std::cerr << "Payload";
    } else if (std::is_same<PayloadType, gbwtgraph::PayloadXL>::value) {
      std::cerr << "PayloadXL";
    } else {
      std::cerr << "Unknown PayloadType";
    }
    std::cerr << std::endl;
  }

  auto start = gbwt::readTimer();

  if (!distance_index) {
    std::function<PayloadType(const pos_t &)> payload_lambda =
        [](const pos_t &pos) -> PayloadType {
      if constexpr (std::is_same<PayloadType, gbwtgraph::PayloadXL>::value) {
        return {0, 0, 0};
      } else {
        return MIPayload::NO_CODE;
      }
    };
    gbwtgraph::index_haplotypes(gbz->graph, index, payload_lambda);
  } else {
    std::function<PayloadType(const pos_t &)> payload_lambda =
        [&](const pos_t &pos) -> PayloadType {
      PayloadType payload = PayloadType::default_payload();
      {
        std::lock_guard<std::mutex> lock(payload_mutex);
        // If we've already seen this node before, return the saved payload
        if (node_id_to_payload.count(id(pos))) {
          payload = node_id_to_payload[id(pos)];
        }
      }
      if (payload != PayloadType::default_payload()) {
        return payload;
      }

      ZipCode zipcode;
      zipcode.fill_in_zipcode(*distance_index, pos);
      gbwtgraph::Payload from_zc = zipcode.get_payload_from_zip();

      if constexpr (std::is_same<PayloadType, gbwtgraph::PayloadXL>::value) {
        payload = PayloadType::from_payload(from_zc);
      } else {
        payload = from_zc;
      }

      if (payload != PayloadType::default_payload()) {
        {
          std::lock_guard<std::mutex> lock(payload_mutex);
          node_id_to_payload.emplace(id(pos), payload);
        }
        return payload;
      } else if (!zipcode_name.empty()) {
        // If they are being saved, add the zipcode to the oversized zipcode
        // list
        zipcode.fill_in_full_decoder();
        {
          std::lock_guard<std::mutex> lock(zipcode_mutex);
          oversized_zipcodes.emplace_back(zipcode);
          size_t zip_index = oversized_zipcodes.size() - 1;
          payload = {0, zip_index};
          node_id_to_payload.emplace(id(pos), payload);
        }
        return payload;
      } else {
        // If the zipcode is too big and we don't have a file to save it
        {
          std::lock_guard<std::mutex> lock(payload_mutex);
          payload = PayloadType::default_payload();
          node_id_to_payload.emplace(id(pos), payload);
        }
        return payload;
      }
    };
    gbwtgraph::index_haplotypes(gbz->graph, index, payload_lambda);
  }

  // Index statistics
  if (show_progress) {
    std::cerr << index.size() << " keys (" << index.unique_keys() << " unique)"
              << std::endl;
    std::cerr << "Minimizer occurrences: " << index.number_of_values()
              << std::endl;
    std::cerr << "Load factor: " << index.load_factor() << std::endl;
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Construction so far: " << seconds << " seconds" << std::endl;
  }

  // Serialize the index
  save_minimizer(index, output_name);

  // If using it, write the larger zipcodes to a file
  if (!zipcode_name.empty()) {
    std::ofstream zip_out(zipcode_name);
    oversized_zipcodes.serialize(zip_out);
    zip_out.close();
  }

  if (show_progress) {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
    std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage())
              << " GiB" << std::endl;
  }
}

template <typename IndexType, typename PayloadType>
size_t
MinimizerIndexHelper<IndexType, PayloadType>::estimate_hash_table_size() const {
  if (show_progress) {
    std::cerr << "Estimating genome size" << std::endl;
  }

  size_t genome_size = 0;

  if (gbz->graph.get_path_count() > 0) {
    gbz->graph.for_each_path_handle([&](const path_handle_t &path_handle) {
      gbz->graph.for_each_step_in_path(
          path_handle, [&](const step_handle_t &step_handle) {
            handle_t handle = gbz->graph.get_handle_of_step(step_handle);
            genome_size += gbz->graph.get_length(handle);
          });
    });
    if (show_progress) {
      std::cerr << "Estimated size based on reference / generic paths: "
                << genome_size << std::endl;
    }
  }

  if (genome_size == 0) {
    gbz->graph.for_each_handle([&](const handle_t &handle) {
      genome_size += gbz->graph.get_length(handle);
    });
    if (show_progress) {
      std::cerr << "Estimated size based on total sequence length: "
                << genome_size << std::endl;
    }
  }

  // Genome size / 2 should be a reasonably tight upper bound for the number of
  // kmers with any specific base in the middle position
  size_t hash_table_size =
      gbwtgraph::KmerIndex<gbwtgraph::Key64, gbwtgraph::Position>::minimum_size(
          genome_size / 2);

  if (show_progress) {
    std::cerr << "Estimated hash table size: 2^"
              << trailing_zeros(hash_table_size) << std::endl;
  }

  return hash_table_size;
}

template <typename IndexType, typename PayloadType>
size_t
MinimizerIndexHelper<IndexType, PayloadType>::trailing_zeros(size_t value) {
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

} // namespace vg

#endif // VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED