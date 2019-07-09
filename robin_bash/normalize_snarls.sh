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
## in terminal:
# gdb vg
#     run mod -F blah test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.vg

export VG_FULL_TRACEBACK=0
set -e

echo compiling!
. ./source_me.sh && make -j 8
echo running!

## testing vg convert
# TEST_CONVERT_FILE=test/robin_haplotypes/test_convert/test_convert
# TEST_CONVERT_FILE=test/robin_chromosomes/test_convert/test_convert
# vg convert -v $TEST_CONVERT_FILE.vg -A >$TEST_CONVERT_FILE.hg
# vg convert -a $TEST_CONVERT_FILE.hg -V >$TEST_CONVERT_FILE.vg

##running normalize_snarls on a full chromosome.
VG_DIR=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10
TEST_DIR=test/robin_chromosomes/chr10
FILE_BASENAME=hgsvc_chr10_construct

# # Jordan's buggy command example
# vg convert -v /public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps/chr10/hgsvc_chr10_construct.vg -A >test/robin_chromosomes/chr10/hgsvc_chr10_construct_test.hg
# echo halfway_there
# vg convert -a test/robin_chromosomes/chr10/hgsvc_chr10_construct_test.hg -P >test/robin_chromosomes/chr10/hgsvc_chr10_construct_test.pg
# echo all_done

# # To produce .snarls:
# vg snarls $VG_DIR/$FILE_BASENAME.vg >$TEST_DIR/$FILE_BASENAME.snarls 
# echo "SNARLS MADE"
# # To produce .gbwt:
# vg index -G $TEST_DIR/$FILE_BASENAME.gbwt -v $VG_DIR/HGSVC.haps.chr10.vcf.gz $VG_DIR/$FILE_BASENAME.vg
# echo "GBWT MADE"
# Convert .vg to .hg:
# vg convert -v $VG_DIR/$FILE_BASENAME.vg -A >$TEST_DIR/$FILE_BASENAME.hg
# echo "CONVERTED VG TO HG"
# Run normalize algorithm:
vg normalize -g $TEST_DIR/$FILE_BASENAME.gbwt -s $TEST_DIR/$FILE_BASENAME.snarls $TEST_DIR/$FILE_BASENAME.hg >$TEST_DIR/$FILE_BASENAME_normalized.hg
echo "NORMALIZED HG MADE"
# # convert .hg to .vg
# vg convert -a $TEST_DIR/$FILE_BASENAME_normalized.hg -V $TEST_DIR/$FILE_BASENAME_normalized.vg
# echo "CONVERTED BACK TO VG."
# # visualize
# ./bin/vg view -dpn $TEST_DIR/$FILE_BASENAME_normalized.vg| \
# dot -Tsvg -o $TEST_DIR/$FILE_BASENAME_normalized.svg
# # chromium-browser $TEST_DIR/$FILE_BASENAME_normalized.svg

# #hg-oriented commands for working on aligning haplotype in middle of snarl. (snarl nodes 23493-23505).
# vg normalize -n test/robin_haplotypes/threads_in_middle_example/chr10_subgraph_0_new.hg >test/robin_haplotypes/threads_in_middle_example/cleaned_mid_hap_snarl.hg
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

