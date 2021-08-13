#!/usr/bin/env bash
#
BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 13

vg msga -f GRCh38_alts/FASTA/HLA/V-352962.fa -t 1 -k 16 | vg mod -U 10 - | vg mod -c - > hla.vg

#flatten it to one path
printf "gi|568815551:1054737-1055734\t0\t1000\n" > region.bed
vg clip hla.vg -b region.bed > clip_flat.vg
vg validate clip_flat.vg
is "$?" 0 "clipped graph is valid"
for i in `vg view clip_flat.vg | grep '^S' | awk '{print $1}'` ; do vg view clip_flat.vg -n $i | grep ^P | grep "gi|568815551:1054737-1055734"; done | wc -l > step_count
vg view clip_flat.vg | grep ^S | wc -l > node_count
diff step_count node_count
is "$?" 0 "every step in clipped graph belongs to reference path"
is $(vg paths -Ev hla.vg -Q "gi|568815551:1054737-1055734" | awk '{ print $2 }') $(vg stats -l clip_flat.vg | awk '{ print $2 }') "clipped graph has same length as ref path"

rm -f region.bed step_count node_count

#clip out one snarl
printf "gi|157734152:29563108-29564082\t90\t92\n" > region.bed
vg clip hla.vg -b region.bed > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is $(vg view clip.vg | grep ^S | wc -l) "49" "Just one node filtered"

rm -f region.bed clip.vg

# clip out one edge
printf "gi|568815564:1054403-1055400\t150\t153\n" > region.bed
vg clip hla.vg -b region.bed > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is $(vg view clip.vg | grep ^L | wc -l) "65" "Just one edge filtered"

rm -f region.bed clip.vg

# clip out low coverage node
vg clip hla.vg -d 4 -P "gi|568815551:1054737-1055734" > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is $(vg view clip.vg | grep ^S | wc -l) "49" "Just one node filtered"

rm -f clip.vg

# clip out out-of-bounds low coverage node
printf "gi|568815551:1054737-1055734\t5\t25\n" > region.bed
vg clip hla.vg -b region.bed -d 4 > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
vg view hla.vg | sort > hla.gfa
vg view clip.vg | sort > clip.gfa
diff hla.gfa clip.gfa
is "$?" 0 "clipping bad region changes nothing"

rm -f clip.vg hla.gfa clip.gfa

# clip out in-bounds low coverage node
printf "gi|568815551:1054737-1055734\t600\t650\n" > region.bed
vg clip hla.vg -b region.bed -d 4 > clip.vg
vg validate clip.vg
is "$?" 0 "clipped graph is valid"
is $(vg view clip.vg | grep ^S | wc -l) "49" "Just one node filtered"

rm -f region.bed clip.vg

rm -f hla.vg


