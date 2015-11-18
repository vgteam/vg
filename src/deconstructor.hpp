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
#include "Fasta.h"
#include "xg.hpp"
namespace vg{
    using namespace std;
    class Deconstructor{
        public:

            Deconstructor();
            ~Deconstructor();
            void clear();
            void set_xg(xg::XG* xg);
            void enumerate_paths_in_index();
            void set_reference(string ref);
            void set_index(Index* index);

            /**
             * Project a path onto another path,
             * much like a transformation of p2 onto p1 in space.
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
            vcflib::Variant mapping_to_variant(Path variant, Path ref);

            /**
             * Turn a vector of variants into a proper VCF.
             */
            vcflib::VariantCallFile write_variants(string filename, vector<vcflib::Variant> variants);
        private:
            // TODO Should probably be able to handle XG or VG indices
            Index* index;
            FastaReference* reference;
            xg::XG* xg;

    };
}
#endif
