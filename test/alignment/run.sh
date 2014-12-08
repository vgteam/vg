#!/bin/bash

if [ $# -lt 5 ]; then
    echo usage: $0 [ref] [num samples] [coverage] [readlength] [fragment length]
    exit
fi

reference=$1
n=$2 # number of samples
coverage=$3
readlength=$4
mfl=$5

echo
echo generating test data
echo

answers=answers.vcf.gz
seed=$(od -An -N2 -i /dev/urandom)
mutatrix -s 0.05 -i 0.01 -g $seed -S sample -p 2 -n $n $reference | bgzip >$answers
tabix -p vcf $answers
samples=`zcat $answers | grep ^#CHROM | cut -f $n- | sed -e "s/\t/\n/g" `
samples=$(for i in $(seq -w $n); do echo sample$i; done)
echo $samples

# calculate coverage
nbp=$(cat $reference | grep -v ">" | sed "s/\n//g" | wc -c)
nreads=$(($nbp * $coverage / $readlength))

refname=`head -1 $reference | sed -e "s/>//"`
for sample in $samples; do cat $sample:$refname:*.fa >$sample.fa; done

for sample in $samples;
do
    seed=$(od -An -N2 -i /dev/urandom)
    mason illumina -n $readlength -N $nreads -rnp $sample -s $seed -sq -mp -ll $mfl -hs 0 -hi 0 $sample.fa
done

