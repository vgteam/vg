#!/bin/bash
date
for burn in 500 999999999
do
    for i in {1..1}
    do 
    (vg mcmc -i 1000 -b $burn -g 100 --vcf-out ./my-output-small/x_small.vcf my-output-small/x_small.mgam ./my-output-small/x_small.vg ./my-output-small/x_small.snarls > ./my-output-small/x_small_paths.vg) 2>&1 | cut -f 1 >> my-output-small/likelihood_${i}_burn${burn}.csv 
    done
done
date