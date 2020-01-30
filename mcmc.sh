set -ex

# Index a couple of nearly identical contigs
# ./bin/vg construct -m 1000 -a -r test/small/x.fa -v test/small/x.vcf.gz >x.vg
# ./bin/vg  index -x x.xg -g x.gcsa -v test/small/x.vcf.gz --gbwt-name x.gbwt -k 16 x.vg
./bin/vg construct -m 1000 -a -r test/1mb1kgp/z.fa -v test/1mb1kgp/z.vcf.gz >x.vg
./bin/vg  index -x x.xg -g x.gcsa -v test/1mb1kgp/z.vcf.gz --gbwt-name x.gbwt -k 16 x.vg

# simulate reads
./bin/vg sim -l 100 -n 100 -e 0.1 -i 0.1 -a -x x.xg > reads.gam

# make the snarls file
./bin/vg snarls x.vg > sites.snarls

# This read is part ref and part alt which matches a haplotype on X, but is possible on Y as well.
./bin/vg mpmap -A -x x.xg -g x.gcsa -G reads.gam > reads.mgam

# call mcmc subcommand
./bin/vg mcmc reads.mgam x.vg sites.snarls > graph_with_paths.vg

# call view on output
# ./bin/vg view -d -w graph_with_paths.vg | dot -Tpdf -o /home/susanna/Desktop/img.pdf
