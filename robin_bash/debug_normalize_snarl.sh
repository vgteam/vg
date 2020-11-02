#!/bin/bash

export VG_FULL_TRACEBACK=1
set -e

echo compiling!
. ./source_me.sh && make -j 8
echo running!

###After Normalization
dir=/home/robin/paten_lab/vg/test/robin_tests/chr21
# base=hgsvc_construct.chr21.robin_made
base=hgsvc_construct.chr21.robin_made.normalized
haps=HGSVC.haps.vcf.gz
threads=16
# meta=200bp.60000num
reads=hgsvc_construct.chr21.robin_made.normalized.read_sim.200bp.60000num.txt
# reads=hgsvc_construct.chr21.robin_made.normalized.read_sim.400bp.5000num.txt

vg map -t $threads -d $dir/$base -x $dir/$base.xg -T $dir/$reads >$dir/$base.alignment.gam
#  vg view -a $dir/$base.alignment.gam -j >$dir/$base.alignment.$meta.json

echo finished mapping.

#jq arbitrary queries for gam files (vg view is easier, more weighty), vg gamcompare with an empty file lets you look at full gam.

###Before Normalization (orignally used in sh file on courtyard - see reconstruct_jmonlong_chr21):
# in_dir=/public/groups/cgl/graph-genomes/jmonlong/hgsvc/haps
# ref=hg38.fa
# in_dir=/home/robin/paten_lab/vg/test/robin_tests/chr21
# vars=HGSVC.haps.vcf.gz
# base=hgsvc_construct.chr21.robin_made
# base_out=hgsvc_construct.chr21.robin_made.test_out
# chrom=chr21
# threads=8

#make graph
#vg construct -r $in_dir/$ref -v $in_dir/$vars -R $chrom -C -m 32 -a -f > $base.vg

#index graph. Note: currently doesn't make .gcsa for some reason.
# vg index -t $threads -x $in_dir/$base_out.xg -G $in_dir/$base_out.gbwt -v $in_dir/$vars -g $in_dir/$base_out.gcsa $in_dir/$base.vg
# vg index -G hgsvc_construct.chr21.robin_made.normalized.gbwt -g hgsvc_construct.chr21.robin_made.normalized.gcsa -v HGSVC.haps.vcf.gz hgsvc_construct.chr21.robin_made.normalized.vg 

#make snarls.
# vg snarls -v $in_dir/$vars $base.vg > $base.snarls.pb

# chunk graph?
# vg chunk -x hgsvc_construct.chr21.robin_made.normalized.xg -G hgsvc_construct.chr21.robin_made.normalized.gbwt -r 0:16608 >chunk_normalized_0_to_16608.vg



###Before and During Normalization
##running normalize_snarls on a full chromosome - local machine.
# TEST_DIR=test/robin_tests/full_chr10
# FILE_NAME=hgsvc_chr10_construct
# FILE_NAME_OUT=junk
# FILE_NAME_OUT=chr10_no_gbwt_handles_at_25128
# FILE_NAME_OUT=hgsvc_chr10_construct_normalized_no_max_size

# TEST_DIR=test/robin_tests/chr21
# FILE_NAME=hgsvc_construct.chr21.robin_made.normalized
# FILE_NAME_OUT=hgsvc_construct.chr21.robin_made.normalized
# FILE_NAME_OUT=hgsvc_construct.chr21.robin_made.normalized.subgraph.301929.exp_context

## running full chr21:
# vg normalize -e -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls.pb $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg

## running subset of chr21
# vg find -x $TEST_DIR/$FILE_NAME.xg -n 301929 -c 100 >$FILE_NAME_OUT.vg
# ./bin/vg view -dpn $FILE_NAME_OUT.vg| \
# dot -Tsvg -o $FILE_NAME_OUT.svg
# chromium-browser $FILE_NAME_OUT.svg


# vg view -dpn $FILE_NAME_OUT.vg| \
# dot -Tsvg -o $FILE_NAME_OUT.svg
# chromium-browser $FILE_NAME_OUT.svg

## for extracting a prenormalized subgraph for looking at chr10

# TEST_DIR=test/robin_tests/chr21
# FILE_NAME=hgsvc_construct.chr21.robin_made
# FILE_NAME_OUT=hgsvc_construct.chr21.robin_made.subgraph.301929.exp_context

# vg find -x $TEST_DIR/$FILE_NAME.xg -n 301929 -c 100 >$FILE_NAME_OUT.vg
# ./bin/vg view -dpn $FILE_NAME_OUT.vg| \
# dot -Tsvg -o $FILE_NAME_OUT.svg
# chromium-browser $FILE_NAME_OUT.svg


##running full chr10
# echo "running normalize (w/ evaluation)"
# valgrind --leak-check=full vg normalize -e -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg
# vg normalize -e -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls.pb $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg


# ##running full chr10 with no max size.
# echo "running normalize (w/ evaluation)"
# vg normalize -e -m 0 -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg


## for printing out the normalized subsnarl:
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >graph_out.vg
# ./bin/vg view -dpn graph_out.vg| \
# dot -Tsvg -o graph_out.svg
# chromium-browser graph_out.svg

## for extracting a prenormalized subgraph for looking at chr10
# vg find -x $TEST_DIR/$FILE_NAME.xg -n 25128 -c 25 >$FILE_NAME_OUT.vg
# ./bin/vg view -dpn $FILE_NAME_OUT.vg| \disambiguating snarl #85 source: 23053 sink: 23075
# dot -Tsvg -o $FILE_NAME_OUT.svg
# chromium-browser $FILE_NAME_OUT.svg

## looking at an old example
# TEST_DIR=test/robin_tests/robin_haplotypes/complex
# FILE_NAME=chr10_subgraph_0_new
# FILE_NAME_OUT=chr10_subgraph_0_new_normalized_200_max_thread_size
# # vg index -G $TEST_DIR/$FILE_NAME.gbwt -v $TEST_DIR/../../HGSVC.haps.chr10.vcf.gz $TEST_DIR/$FILE_NAME.vg
# # vg convert -v $TEST_DIR/$FILE_NAME.vg -A >$TEST_DIR/$FILE_NAME.hg
# vg normalize -e -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg
# vg convert -a $TEST_DIR/$FILE_NAME_OUT.hg -V >$TEST_DIR/$FILE_NAME_OUT.vg
# vg mod -g 609548 -x 65 $TEST_DIR/$FILE_NAME_OUT.vg | vg view -dpn - | dot -Tsvg -o $TEST_DIR_OUT/$FILE_NAME_OUT.svg
# chromium-browser $TEST_DIR_OUT/$FILE_NAME_OUT.svg




### After Normalization:
## for making a snarls file:
# vg convert -a $TEST_DIR/$FILE_NAME_OUT.hg -V >$TEST_DIR/$FILE_NAME_OUT.vg
# echo "hg converted to vg"
# vg snarls $TEST_DIR/$FILE_NAME_OUT.vg >$TEST_DIR/$FILE_NAME_OUT.snarls 
# echo ".snarls made"

## for evaluating normalized graph:
# echo "getting vg stats:"
# vg stats -z -l $TEST_DIR/$FILE_NAME_OUT.vg

## creating a new gbwt graph from the outgraph:
# vg index -G $TEST_DIR/$FILE_NAME_OUT.gbwt -v $TEST_DIR/../HGSVC.haps.chr10.vcf.gz $TEST_DIR/$FILE_NAME_OUT.vg
# echo "gbwt made"

