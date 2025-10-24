/*
minimizer_index_helper.hpp
Shared logic for minimizer index construction.
*/

#pragma once
#include <fstream>
#include <functional>
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

template <typename IndexType, typename PayloadType>
void build_minimizer_index(const gbwtgraph::GBZ *gbz, IndexType &index,
                           const bdsg::SnarlDistanceIndex &distance_index,
                           const std::string &distance_name,
                           const std::string &zipcode_name,
                           const std::string &output_name,
                           size_t hash_table_size, size_t threshold,
                           size_t iterations, bool space_efficient_counting,
                           bool weighted, bool use_syncmers,
                           bool use_distance_index, bool use_zipcode_index,
                           bool progress) {

  // Zipcodes
  // oversized_zipcodes may be stored alongside the minimizer index in the file
  // specified by zipcode_name
  ZipCodeCollection oversized_zipcodes;

  // Map node id to what gets stored in the payload - either the zipcode or
  // index into oversized_zipcodes
  hash_map<vg::id_t, PayloadType> node_id_to_payload;
  node_id_to_payload.reserve(gbz->graph.max_node_id() -
                             gbz->graph.min_node_id());

  // Build the index.
  if (progress) {
    std::cerr << "Building MinimizerIndex with k = " << index.k();
    if (index.uses_syncmers()) {
      std::cerr << ", s = " << index.s();
    } else {
      std::cerr << ", w = " << index.w();
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
  if (!use_distance_index) {
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
#pragma omp critical
      {
        // If we've already seen this node before, then return the saved payload
        if (node_id_to_payload.count(id(pos))) {
          payload = node_id_to_payload[id(pos)];
        }
      }
      if (payload != PayloadType::default_payload()) {
        return payload;
      }

      ZipCode zipcode;
      zipcode.fill_in_zipcode(distance_index, pos);
      gbwtgraph::Payload from_zc = zipcode.get_payload_from_zip();

      if constexpr (std::is_same<PayloadType, gbwtgraph::PayloadXL>::value) {
        payload = PayloadType::from_payload(from_zc);
      } else {
        payload = from_zc;
      }

      if (payload != PayloadType::default_payload()) {
#pragma omp critical
        {
          node_id_to_payload.emplace(id(pos), payload);
        }
        return payload;
      } else if (use_zipcode_index) {
        // Otherwise, if they are being saved, add the zipcode to the oversized
        // zipcode list And remember the zipcode

        // Fill in the decoder to be saved too
        zipcode.fill_in_full_decoder();

#pragma omp critical
        {
          oversized_zipcodes.emplace_back(zipcode);
          size_t zip_index = oversized_zipcodes.size() - 1;
          payload = {0, zip_index};
          node_id_to_payload.emplace(id(pos), payload);
        }
        return payload;
      } else {
// If the zipcode is too big and we don't have a file to save the big zipcodes
#pragma omp critical
        {
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
  if (use_zipcode_index) {
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