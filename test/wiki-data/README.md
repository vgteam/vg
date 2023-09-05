# Test Data Directory
This directory contains mock data for testing the examples on the VG wiki.
It has the same shape as some whole-genome data, but is not actually whole-genome scale.

# mock-hs37d5

Can be generated like this:

        
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa ~/build/vg/trash/hs37d5.fa
    grep -v "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" ~/build/vg/trash/hs37d5.fa | grep "^>" -A20 | grep -v "^--" > ~/build/vg/trash/hs37d5.part.fa
   
   
    mv hs37d5.fa hs37d5.orig.fa
    mv hs37d5.part.fa hs37d5.fa
    

    for URL in  https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz ; do
        BASENAME=$(basename ${URL})
        # Expect Curl to complain when head stops reading
        curl -sSL ${URL} | zcat - | head -n2000 | grep "^#" > ${BASENAME%.gz}
    done
    
    for CHR in {1..22} ; do
        mv ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf chr${CHR}.vcf
    done
    
    mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf chrX.vcf
    mv ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf chrY.vcf
    mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf chrMT.vcf
                  
    for VCF in *.vcf ; do
        bgzip ${VCF}
        tabix -p vcf ${VCF}.gz
    done
    
Simulated reads can then be generated like this:

    vg construct -r hs37d5.fa >graph.vg
    vg sim -x graph.vg -l 150 -p 500 -v 5 -n 10000 -P MT -a >sim.gam
    vg view -aXi sim.gam >sim.fq
    
Since the contigs are short, we probably won't actually get 500 bp read separation.
