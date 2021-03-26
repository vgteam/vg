#!/bin/bash
## to give permission to run script:
# chmod +x test_F_on_x.sh
## this bash to be run on dir /vg/

# VG=chr10_subgraph_2dels-shift-729006
# VG=clean_chr10_subgraph_0_new
# VG=chr10_subgraph_2dels-323159
# VG=chr10_subgraph_3ins-1558671

## useful debug tools:
# export VG_FULL_TRACEBACK=1
    
# valgrind vg mod -F blah test/robin_snarl_examples/chr10_subgraph_0_new.vg 
    #Note: to make more informative, commment out the two lines in MakeFile under "use 
    #jemalloc", delete bin/vg, and recompile.

## in terminal:
# gdb vg
#     run normalize -g test/robin_tests/full_chr10/hgsvc_chr10_construct.gbwt -s test/robin_tests/full_chr10/hgsvc_chr10_construct.snarls test/robin_tests/full_chr10/hgsvc_chr10_construct.hg >test/robin_tests/full_chr10/hgsvc_chr10_construct_normalized.hg

export VG_FULL_TRACEBACK=1
set -e

echo compiling!
. ./source_me.sh && make -j 8
echo running!

## constructing a smaller graph from a larger one - one method of subsetting a graph.
# vg construct -r test/small/x.fa -v test/small/x.vcf.gz -R x:1-10

# TEST_DIR=test/robin_tests/vis_vg_find_sample
# FILE_NAME=chr10_subset_vg_find

# ./bin/vg view -dpn $TEST_DIR/$FILE_NAME.vg| \
# dot -Tsvg -o $TEST_DIR/$FILE_NAME.svg
# chromium-browser $TEST_DIR/$FILE_NAME.svg

##running normalize_snarls on a full chromosome - local machine.
TEST_DIR=test/robin_tests/full_chr10
FILE_NAME=hgsvc_chr10_construct

## for printing out the subsnarl:
vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/full_chr10_normalized.vg
# visualize
# ./bin/vg view -dpn $TEST_DIR/graph3.vg| \
# dot -Tsvg -o $TEST_DIR/graph3.svg
# chromium-browser $TEST_DIR/graph3.svg

# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/hgsvc_chr10_construct_normalized_one_snarl_pre_long_chr10.hg
# echo "NORMALIZED HG MADE"
# # convert .hg to .vg
# vg convert -a $TEST_DIR/hgsvc_chr10_construct_normalized_one_snarl_pre_long_chr10.hg -V $TEST_DIR/hgsvc_chr10_construct_normalized_one_snarl_pre_long_chr10.vg
# echo "CONVERTED BACK TO VG."
# #extract subsnarl:
# vg find 

# # visualize
# ./bin/vg view -dpn hgsvc_chr10_one_snarl_pre_long_chr10_extracted.vg| \
# dot -Tsvg -o hgsvc_chr10_one_snarl_pre_long_chr10_extracted.svg
# chromium-browser hgsvc_chr10_one_snarl_pre_long_chr10_extracted.svg






# echo "CONVERTED VG TO HG"
# Run normalize algorithm:
# valgrind --leak-check=full vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/hgsvc_chr10_construct_normalized.hg
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >big_chr10_path.vg
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/hgsvc_chr10_construct_normalized.hg
# echo "NORMALIZED HG MADE"
# # convert .hg to .vg
# vg convert -a $TEST_DIR/hgsvc_chr10_construct_normalized.hg -V $TEST_DIR/hgsvc_chr10_construct_normalized.vg
# echo "CONVERTED BACK TO VG."
# # visualize
# ./bin/vg view -dpn $TEST_DIR/$FILE_NAME.vg| \
# dot -Tsvg -o $TEST_DIR/$FILE_NAME.svg
# chromium-browser $TEST_DIR/$FILE_NAME.svg
# ./bin/vg view -dpn big_chr10_path.vg| \
# dot -Tsvg -o big_chr10_path.svg
# chromium-browser big_chr10_path.svg


# # visualize subsetted chr10
# vg mod -g 7280 -x 5360 $TEST_DIR/$FILE_NAME.vg >$TEST_DIR/hgsvc_chr10_construct_first_few_snarls.vg #produces snarl from 1879:12785
# echo "subgraph made"
# # vg find -x $TEST_DIR/hgsvc_chr10_construct.xg -n 7280 -c 5360 >$TEST_DIR/$FILE_NAME.vg #1878:12785 node range.
# ./bin/vg view -dpn $TEST_DIR/hgsvc_chr10_construct_first_few_snarls.vg| \
# dot -Tsvg -o $TEST_DIR/hgsvc_chr10_construct_first_few_snarls.svg
# chromium-browser $TEST_DIR/$FILE_BASENAME_normalized.svg

# ## split off the first few snarls from chromosome ten: (aiming for nodes between 1883 and 12677)
# # VG_DIR=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10
# TEST_DIR=test/robin_tests/chr10_subset/set_1
# FILE_NAME=chr10_subset_vg_find
# # vg find -x $VG_DIR/hgsvc_chr10_construct.xg -p "chr10:1883-12677" -c 10 >$TEST_DIR/$FILE_NAME.vg
# vg find -x $TEST_DIR/hgsvc_chr10_construct.xg -n 7280 -c 5360 >$TEST_DIR/$FILE_NAME.vg #1878:12785 node range.
# # vg mod -g 3000 -x 5 $VG_DIR/hgsvc_chr10_construct.vg >$TEST_DIR/$FILE_NAME.vg 
# echo "vg subgraph made!"

# ## split off the first few snarls from chromosome ten: (aiming for nodes between 1883 and 12677)
# VG_DIR=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10
# TEST_DIR=test/robin_tests/chr10_subset/set_1
# FILE_NAME=chr10_subset_vg_find
# # vg find -x $VG_DIR/hgsvc_chr10_construct.xg -p "chr10:1883-12677" -c 10 >$TEST_DIR/$FILE_NAME.vg
# vg find -x $VG_DIR/hgsvc_chr10_construct.xg -n 7280 -c 5360 >$TEST_DIR/$FILE_NAME.vg #1878:12785 node range.
# # vg mod -g 3000 -x 5 $VG_DIR/hgsvc_chr10_construct.vg >$TEST_DIR/$FILE_NAME.vg 
# echo "vg subgraph made!"


# ## run normalize_snarls on subsetted full chromosome 10:
# TEST_DIR=test/robin_tests/chr10_subset/set_1
# FILE_NAME=chr10_subset_vg_find
# # To produce .snarls:
# vg snarls $TEST_DIR/$FILE_NAME.vg >$TEST_DIR/$FILE_NAME.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/$FILE_NAME.gbwt -v $TEST_DIR/../../HGSVC.haps.chr10.vcf.gz $TEST_DIR/$FILE_NAME.vg
# echo "GBWT MADE"
# # Convert .vg to .hg:
# vg convert -v $TEST_DIR/$FILE_NAME.vg -A >$TEST_DIR/$FILE_NAME.hg
# echo "CONVERTED VG TO HG"
# # vg convert -a $TEST_DIR/$FILE_NAME.hg -V >$TEST_DIR/$FILE_NAME.vg
# # echo "CONVERTED HG TO VG"
# # Run normalize algorithm:
# ls $TEST_DIR
# # echo $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg $TEST_DIR/chr10_subset_normalized.hg
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/chr10_subset_vg_find_normalized.hg
# echo "normalized."


# ## run normalize_snarls on local machine small test:
# TEST_DIR=test/robin_tests/test_normalize
# FILE_NAME=chr10_subgraph_0_new
# # To produce .snarls:
# vg snarls $TEST_DIR/$FILE_NAME.vg >$TEST_DIR/$FILE_NAME.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/$FILE_NAME.gbwt -v $TEST_DIR/../HGSVC.haps.chr10.vcf.gz $TEST_DIR/$FILE_NAME.vg
# echo "GBWT MADE"
# # Convert .vg to .hg:
# vg convert -v $TEST_DIR/$FILE_NAME.vg -A >$TEST_DIR/$FILE_NAME.hg
# # vg convert -v $TEST_DIR/chr10_subgraph_0_new.vg -A >$TEST_DIR/$FILE_NAME.hg
# echo "CONVERTED VG TO HG"
# # Run normalize algorithm:
# echo $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg $TEST_DIR/chr10_subgraph_0_new_normalized.hg
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/chr10_subgraph_0_new_normalized.hg
# echo "normalized."

# ## split off the first few snarls from chromosome ten: (aiming for nodes between 1883 and 12677)
# VG_DIR=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10
# TEST_DIR=test/robin_tests/chr10_subset
# vg mod -g 7280 -x 5360 $VG_DIR/hgsvc_chr10_construct.vg >$TEST_DIR/hgsvc_chr10_construct_first_few_snarls.vg #produces snarl from 1879:12785

# ## run normalize_snarls on already made full chromosome 10:
# TEST_DIR=test/robin_tests/chr10
# vg normalize -g $TEST_DIR/hgsvc_chr10_construct.gbwt -s $TEST_DIR/hgsvc_chr10_construct.snarls $TEST_DIR/hgsvc_chr10_construct.hg >$TEST_DIR/$FILE_BASENAME_normalized.hg
# echo "NORMALIZED HG MADE"

# ##running normalize_snarls on a full chromosome.
# VG_DIR=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10
# TEST_DIR=test/robin_tests/chr10
# FILE_BASENAME=hgsvc_chr10_construct
# # To produce .snarls:
# vg snarls $VG_DIR/$FILE_BASENAME.vg >$TEST_DIR/$FILE_BASENAME.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/$FILE_BASENAME.gbwt -v $VG_DIR/HGSVC.haps.chr10.vcf.gz $VG_DIR/$FILE_BASENAME.vg
# echo "GBWT MADE"
# # Convert .vg to .hg:
# vg convert -v $VG_DIR/$FILE_BASENAME.vg -A >$TEST_DIR/$FILE_BASENAME.hg
# echo "CONVERTED VG TO HG"
# Run normalize algorithm:
# vg normalize -g $TEST_DIR/$FILE_BASENAME.gbwt -s $TEST_DIR/$FILE_BASENAME.snarls $TEST_DIR/$FILE_BASENAME.hg >$TEST_DIR/$FILE_BASENAME_normalized.hg
# echo "NORMALIZED HG MADE"
# # convert .hg to .vg
# vg convert -a $TEST_DIR/$FILE_BASENAME_normalized.hg -V $TEST_DIR/$FILE_BASENAME_normalized.vg
# echo "CONVERTED BACK TO VG."
# # visualize
# ./bin/vg view -dpn $TEST_DIR/$FILE_BASENAME_normalized.vg| \
# dot -Tsvg -o $TEST_DIR/$FILE_BASENAME_normalized.svg
# # chromium-browser $TEST_DIR/$FILE_BASENAME_normalized.svg


# ## testing vg normalize in Courtyard on smaller graph:
# TEST_DIR=test/robin_tests/normalize_2
# # To produce .snarls:
# vg snarls $TEST_DIR/normalize.vg >$TEST_DIR/normalize.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/normalize.gbwt -v $TEST_DIR/HGSVC.haps.chr10.vcf.gz $TEST_DIR/normalize.vg
# echo "GBWT MADE"
# # Convert .vg to .hg:
# vg convert -v $TEST_DIR/normalize.vg -A >$TEST_DIR/normalize.hg
# echo "CONVERTED VG TO HG"
# # Run normalize algorithm:
# vg normalize -g $TEST_DIR/normalize.gbwt -s $TEST_DIR/normalize.snarls $TEST_DIR/normalize.hg >$TEST_DIR/normalize_out.hg
# echo "normalized."

# ## testing vg normalize in local machine on smaller graph (checking that serialization still works):
# TEST_DIR=test/robin_tests/robin_haplotypes/threads_in_middle_example
# # To produce .snarls:
# vg snarls $TEST_DIR/chr10_subgraph_0_new.vg >$TEST_DIR/new_remake_test/normalize.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/new_remake_test/normalize.gbwt -v $TEST_DIR/HGSVC.haps.chr10.vcf.gz $TEST_DIR/chr10_subgraph_0_new.vg
# echo "GBWT MADE"
# # Convert .vg to .hg:
# vg convert -v $TEST_DIR/chr10_subgraph_0_new.vg -A >$TEST_DIR/new_remake_test/normalize.hg
# # vg convert -v $TEST_DIR/chr10_subgraph_0_new.vg -A >$TEST_DIR/new_remake_test/normalize.hg
# echo "CONVERTED VG TO HG"
# # Run normalize algorithm:
# vg normalize -g $TEST_DIR/new_remake_test/normalize.gbwt -s $TEST_DIR/new_remake_test/normalize.snarls $TEST_DIR/new_remake_test/normalize.hg >$TEST_DIR/new_remake_test/normalize_out.hg
# echo "normalized."

# ## Jordan's buggy command example
# vg convert -v /public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10/hgsvc_chr10_construct.vg -A >test/robin_tests/chr10/hgsvc_chr10_construct_test.hg
# echo halfway_there
# vg convert -a test/robin_tests/chr10/hgsvc_chr10_construct_test.hg -P >test/robin_tests/chr10/hgsvc_chr10_construct_test.pg
# echo all_done
# ## testing vg convert
# TEST_CONVERT_FILE=test/robin_haplotypes/test_convert/test_convert
# TEST_CONVERT_FILE=test/robin_chromosomes/test_convert/test_convert
# vg convert -v $TEST_CONVERT_FILE.vg -A >$TEST_CONVERT_FILE.hg
# vg convert -a $TEST_CONVERT_FILE.hg -V >$TEST_CONVERT_FILE.vg

# ##hg-oriented commands for working on aligning haplotype in middle of snarl. (snarl nodes 23493-23505).
# TEST=test/robin_tests/robin_haplotypes/threads_in_middle_example
# vg normalize -g $TEST/chr10_subgraph_0_new.gbwt -s $TEST/chr10_subgraph_0_new.snarls $TEST/chr10_subgraph_0_new.hg >$TEST/cleaned_mid_hap_snarl.hg
# vg convert -a test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.hg -V >test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl_from_hash.vg
# ./bin/vg view -dpn test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl_from_hash.vg |
#     dot -Tsvg -o test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl_from_hash.svg
# chromium-browser test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl_from_hash.svg

# #for working on aligning haplotype in middle of snarl. (snarl nodes 23493-23505). [the last use of mod_main before using normalize_main].
# ./bin/vg mod -F blah test/robin_haplotypes/threads_in_middle_example/chr10_subgraph_0_new.vg >test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.vg
# ./bin/vg view -dpn test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.vg| \
# dot -Tsvg -o test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.svg
# chromium-browser test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.svg

# vg mod -F blah test/robin_snarl_examples/chr10_subgraph_0_new.vg

# To produce a gbwt file:
# vg index -x chr10_subgraph_0_new.xg -G chr10_subgraph_0_new.gbwt -v HGSVC.haps.chr10.vcf.gz chr10_subgraph_0_new.vg

# vg mod -F blah test/robin_snarl_examples/chr10_subgraph_0_new.vg >test/robin_snarl_examples/cleaned_snarl.vg
# # valgrind vg mod -F blah test/robin_snarl_examples/chr10_subgraph_0_new.vg
# bin/vg view -dpn test/robin_snarl_examples/chr10_subgraph_0_new.vg| \
# dot -Tsvg -o test/robin_snarl_examples/chr10_subgraph_0_new_1.svg
# chromium-browser test/robin_snarl_examples/chr10_subgraph_0_new_1.svg

# vg mod -F blah test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.vg
# vg mod -F blah test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.vg >test/robin_haplotypes/simple/modified_snarl/modified_snarl.vg
# bin/vg view -dpn test/robin_haplotypes/simple/modified_snarl/modified_snarl.vg| \
# dot -Tsvg -o test/robin_haplotypes/simple/modified_snarl/modified_snarl.svg
# chromium-browser test/robin_haplotypes/simple/modified_snarl/modified_snarl.svg

# gdb vg mod -F blah test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.vg

##for looking at snarl number:
# . ./source_me.sh && make -j 8
# bin/vg snarls test/$VG.vg > test/$VG.snarls
# bin/vg mod -F test/$VG.snarls test/$VG.vg

## for running clean_all_snarls:
# . ./source_me.sh && make -j 8
# bin/vg snarls test/$VG.vg > test/$VG.snarls
# bin/vg mod -F test/$VG.snarls test/$VG.vg >test/clean_$VG.vg
# bin/vg view -dpn test/clean_$VG.vg| \
# dot -Tsvg -o test/clean_$VG.svg
# chromium-browser test/clean_$VG.svg

# . ./source_me.sh && make -j 8
# bin/vg mod -F "220 218" test/0_cluttered_snarl_simple.vg # for misc. terminal output (0_demo_final_0)
# bin/vg mod -F "1 6" test/bash_x_out_1.vg # for misc. terminal output (old)

##View "cluttered" snarl for 0_demo_final_0:
# bin/vg view -dpn test/0_cluttered_snarl_simple.vg| \
# dot -Tsvg -o test/0_cluttered_snarl_simple.svg
# chromium-browser test/0_cluttered_snarl_simple.svg

##Testing 0_demo_final_0 on 0_cluttered_snarl_simple.vg:
# bin/vg mod -F "220 218" test/0_cluttered_snarl_simple.vg >test/0_clean_snarl_simple.vg
# bin/vg view -dpn test/0_clean_snarl_simple.vg| \
# dot -Tsvg -o test/0_clean_snarl_simple.svg
# chromium-browser test/0_clean_snarl_simple.svg

##Testing 0_demo_final_0 clean_all_snarls on 0_cluttered_snarl_simple.vg with snarl:
# bin/vg snarls test/0_cluttered_snarl_simple.vg > test/0_cluttered_snarl_simple.snarls

##view terminal output:
# bin/vg mod -F "test/0_cluttered_snarl_simple.snarls" test/0_cluttered_snarl_simple.vg

##view graph:
# bin/vg mod -F "test/0_cluttered_snarl_simple.snarls" test/0_cluttered_snarl_simple.vg >test/0_clean_all_snarl_simple.vg
# bin/vg view -dpn test/0_clean_all_snarl_simple.vg| \
# dot -Tsvg -o test/0_clean_all_snarl_simple.svg
# chromium-browser test/0_clean_all_snarl_simple.svg

##view snarls:
# vg view -Rj test/0_cluttered_snarl_simple.snarls

## The following was for when F was adding 6 T's to all the nodes in the first snarl.
## its output was bash_x_out_1.vg and bash_x_out_1.svg:

# bin/vg mod -F "1 6" test/x.vg >test/bash_x_out_1.vg
# bin/vg view -dpn test/bash_x_out_1.vg| \
# dot -Tsvg -o test/bash_x_out_1.svg

# chromium-browser test/bash_x_out_1.svg

## Now, using graph bash_x_out.vg, consolidate duplicated T's into a single node.

# bin/vg mod -F "1 6" test/bash_x_out_1.vg >test/bash_x_out_clean_snarl.vg
# bin/vg view -dpn test/bash_x_out_1_clean_snarl.vg| \
# dot -Tsvg -o test/bash_x_out_1_clean_snarl.svg

# chromium-browser test/bash_x_out_1_clean_snarl.svg

## testing 0_demo_align_strings.cpp

# bin/vg mod -F "1 6" test/bash_x_out_1.vg >test/0_align_strings.vg
# bin/vg view -dpn test/0_align_strings.vg| \
# dot -Tsvg -o test/0_align_strings.svg

# chromium-browser test/0_align_strings.svg

# bin/vg mod -F "1 6" test/bash_x_out_1.vg >test/0_align_strings.vg
# bin/vg view -dpn test/0_align_strings.vg| \
# dot -Tsvg -o test/0_align_strings.svg

# chromium-browser test/0_align_strings.svg

# vg view -dpn chr10_subgraph_0_new.vg| \
# dot -Tsvg -o chr10_subgraph_0_new_2.svg
# chromium-browser chr10_subgraph_0_new_2.svg

