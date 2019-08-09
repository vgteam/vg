
export VG_FULL_TRACEBACK=1
set -e

echo compiling!
. ./source_me.sh && make -j 8
echo running!

###Before and During Normalization
##running normalize_snarls on a full chromosome - local machine.
TEST_DIR=test/robin_tests/full_chr10
FILE_NAME=hgsvc_chr10_construct
FILE_NAME_OUT=hgsvc_chr10_construct_normalized_no_max_size
# FILE_NAME_OUT=junk

##running full chr10
echo "running normalize (w/ evaluation)"
vg normalize -e -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >$TEST_DIR/$FILE_NAME_OUT.hg

## for printing out the normalized subsnarl:
# vg normalize -g $TEST_DIR/$FILE_NAME.gbwt -s $TEST_DIR/$FILE_NAME.snarls $TEST_DIR/$FILE_NAME.hg >graph_out.vg
# ./bin/vg view -dpn graph_out.vg| \
# dot -Tsvg -o graph_out.svg
# chromium-browser graph_out.svg

## for extracting a prenormalized subgraph for looking at chr10
# vg find -x test/robin_tests/full_chr10/hgsvc_chr10_construct.xg -n 1236806 -c 25 >graph_in.vg
# ./bin/vg view -dpn graph_in.vg| \
# dot -Tsvg -o graph_in.svg
# chromium-browser graph_in.svg

### After Normalization:
## for making a snarls file:
vg convert -a $TEST_DIR/$FILE_NAME_OUT.hg -V >$TEST_DIR/$FILE_NAME_OUT.vg
echo "hg converted to vg"
vg snarls $TEST_DIR/$FILE_NAME_OUT.vg >$TEST_DIR/$FILE_NAME_OUT.snarls 
echo ".snarls made"

## for evaluating normalized graph:
echo "getting vg stats:"
vg stats -z -l $TEST_DIR/$FILE_NAME_OUT.vg

## creating a new gbwt graph from the outgraph:
vg index -G $TEST_DIR/$FILE_NAME_OUT.gbwt -v $TEST_DIR/../HGSVC.haps.chr10.vcf.gz $TEST_DIR/$FILE_NAME_OUT.vg
echo "gbwt made"