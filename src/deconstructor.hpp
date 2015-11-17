#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "Variant.h"
#include "index.hpp"
#include "path.hpp"
#include "vg.hpp"
#include "vg.pb.h"
namespace vg{
    using namespace std;
    class Deconstructor{
        public:

            Deconstructor();
            ~Deconstructor();
            void clear();
            void set_index(Index ind, bool isXG);

            /**
             * Project a mapping onto another mapping.
             */
            Path relative_mapping(Path& p1, Path& p2);

            /**
             * Convenience function.
             * Project a mapping with name 'p2' onto p1.
             * Return the projection as a map.
             */
            Path relative_mapping(Path& p1, string p2);

            /**
             * Build a vcf record from two paths, with the
             * second path argument taken as the reference.
             */
            Variant mapping_to_variant(Path variant, Path ref);

            /**
             * Turn a vector of variants into a proper VCF.
             */
            VariantCallFile write_variants(string filename, vector<Variant> variants);
        private:
            // TODO Should probably be able to handle XG or VG indices
            Index index;

    }
}

#endif
