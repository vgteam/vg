set -ex

# Index a couple of nearly identical contigs
vg construct -m 1000 -a -r test/small/x.fa -v test/small/x.vcf.gz >x.vg
vg index -x x.xg -g x.gcsa -v test/small/x.vcf.gz --gbwt-name x.gbwt -k 16 x.vg

# simulate reads
vg sim -l 100 -n 100 -e 0.1 -i 0.1 -a -x x.xg > reads.gam

# make the snarls file
vg snarls x.vg > sites.snarls

# This read is part ref and part alt which matches a haplotype on X, but is possible on Y as well.
vg mpmap -A -x x.xg -g x.gcsa -G reads.gam > reads.mgam

# call mcmc subcommand
vg mcmc reads.mgam x.vg sites.snarls > graphs_with_paths.vg