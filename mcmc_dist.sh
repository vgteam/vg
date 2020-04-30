#!bin/bash
start=`date +%T`

for i in {1..100}
do 
    (vg mcmc -i 10000 --vcf-out ./my-output-small/x_small.vcf my-output-small/x_small.mgam ./my-output-small/x_small.vg ./my-output-small/x_small.snarls > ./my-output-small/x_small_paths.vg) 2>&1 | cut -f 1 >> my-output-small/likelihood_${i}.csv 
    hap.py test/small/x.vcf.gz my-output-small/x_small.vcf -o my-output-small/hapPy_results --force-interactive
    cut -d "," -f 1,13 my-output-small/hapPy_results.summary.csv | uniq | cat > my-output-small/dist.csv
    grep "SNP" my-output-small/dist.csv | cut -d "," -f 2 - >> my-output-small/SNP_dist.csv
    grep "INDEL" my-output-small/dist.csv | cut -d "," -f 2 - >> my-output-small/INDEL_dist.csv 
done
end=`date +%T`

echo $start
echo $end 

#grep "clikelihood" my-output-small/likelihood_6.csv | cut -d " " -f 2 - >> my-output-small/likelihood_extracted_6.csv


