/**
\file minimizer_index_helper.hpp
Shared logic for minimizer index construction.
*/
#ifndef VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED
#define VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED

#include <fstream>
#include <functional>
#include <mutex>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/index.h>
#include <omp.h>
#include <string>
#include <type_traits>
#include <vector>

#include "snarl_distance_index.hpp"
#include "vg/io/vpkg.hpp"
#include "zip_code.hpp"

namespace mi_helper {

using namespace std;
using namespace vg;

/**
 * Count the number of trailing zero bits in a value.

 * @return Number of trailing zeros, or 0 if value is 0
 */
size_t trailing_zeros(size_t value);

/**
 * Estimate an appropriate hash table size for a minimizer index based on genome size.
 * 
 * Calculates genome size from paths if available, otherwise from total sequence length.
 * Returns a power-of-2 size suitable for storing kmers.
 *
 * @return Estimated hash table size (always a power of 2)
 */
size_t estimate_hash_table_size(const gbwtgraph::GBZ &gbz, bool progress);

/**
 * Identify and mark frequent kmers in a minimizer index.
 * 
 * Finds kmers that occur more than the specified threshold and adds them
 * to the index's frequent kmer list for special handling.
 */
template <typename IndexType>
void set_frequent_kmers(const gbwtgraph::GBZ *gbz, IndexType &index, int k,
                        size_t threshold, bool space_efficient_counting,
                        size_t hash_table_size, size_t iterations,
                        bool progress) {

  vector<gbwtgraph::Key64> frequent_kmers;
  if (progress) {
    string algorithm = (space_efficient_counting ? "space-efficient" : "fast");
    cerr << "Finding frequent kmers using the " << algorithm << " algorithm"
         << endl;
  }

  if (hash_table_size == 0) {
    hash_table_size = estimate_hash_table_size(*gbz, progress);
  }
  frequent_kmers = gbwtgraph::frequent_kmers<gbwtgraph::Key64>(
      gbz->graph, k, threshold, space_efficient_counting, hash_table_size);
  if (progress) {
    cerr << "Found " << frequent_kmers.size() << " kmers with more than "
         << threshold << " hits" << endl;
  }
  if (!frequent_kmers.empty()) {
    index.add_frequent_kmers(frequent_kmers, iterations);
  }
}

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
template <typename IndexType, typename PayloadType>
void build_minimizer_index(const gbwtgraph::GBZ *gbz, IndexType &index,
                           const bdsg::SnarlDistanceIndex *distance_index,
                           const std::string &zipcode_name,
                           const std::string &output_name, bool progress) {

  // Zipcodes
  // oversized_zipcodes may be stored alongside the minimizer index in the file
  // specified by zipcode_name
  ZipCodeCollection oversized_zipcodes;

  // Map node id to what gets stored in the payload - either the zipcode or
  // index into oversized_zipcodes
  vg::hash_map<vg::id_t, PayloadType> node_id_to_payload;
  node_id_to_payload.reserve(gbz->graph.max_node_id() -
                             gbz->graph.min_node_id());

  // Mutex for protecting shared data structures in parallel regions
  std::mutex payload_mutex;
  std::mutex zipcode_mutex;

  // Build the index.
  if (progress) {
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
        // If we've already seen this node before, then return the saved payload
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
        // Otherwise, if they are being saved, add the zipcode to the oversized
        // zipcode list And remember the zipcode

        // Fill in the decoder to be saved too
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
        // If the zipcode is too big and we don't have a file to save the big zipcodes
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
  // Index statistics.
  if (progress) {
    std::cerr << index.size() << " keys (" << index.unique_keys() << " unique)"
              << std::endl;
    std::cerr << "Minimizer occurrences: " << index.number_of_values()
              << std::endl;
    std::cerr << "Load factor: " << index.load_factor() << std::endl;
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Construction so far: " << seconds << " seconds" << std::endl;
  }

  // Serialize the index.
  save_minimizer(index, output_name);

  // If using it, write the larger zipcodes to a file
  if (!zipcode_name.empty()) {
    ofstream zip_out(zipcode_name);
    oversized_zipcodes.serialize(zip_out);
    zip_out.close();
  }

  if (progress) {
    double seconds = gbwt::readTimer() - start;
    std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
    std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage())
              << " GiB" << std::endl;
  }
}

} // namespace mi_helper

#endif // VG_MINIMIZER_INDEX_HELPER_HPP_INCLUDED