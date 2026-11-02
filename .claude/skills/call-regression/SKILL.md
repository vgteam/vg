---
name: call-regression
description: Regression-test vg call's read-likelihood caller (PR 4990's --read-likelihood, nesting, phasing, anchors and mosaic) on real HG002 reads, short and ONT, on chr20 or chr6, by checking that the VCF, mosaic and anchors are byte-identical to a baseline and by scoring them against GIAB v5.0q. Use before pushing or merging any change that could touch vg call's output, and to set a new baseline after an intended output change.
---

# Regression-test vg call on HG002 chr20 and chr6

The gate calls HG002 on one contig twice with the vg under test: once from NovaSeq 30x short
reads and once from E821 ONT reads. Each read set has its own sampled graph. The gate then
checks the output files against the baseline's md5s and scores the calls against GIAB's
HG002 v5.0q benchmark.

There are baselines for two contigs. chr20 is the quick check, at about 6 minutes. chr6 is
the larger one, at about 30 minutes including scoring. It has 2.7 times as many records,
and the MHC, which is one of the hardest regions to call. Set `CONTIG=chr6` for every
script to use it; the default is `CONTIG=chr20`.

For a change that should not move the output, such as a refactor, a speed-up or a change to
another subcommand, every file must be byte-identical, and that is the whole test. For a
change that is meant to move the output, the scores say whether accuracy rose or fell. Once
the change is accepted, its run becomes the new baseline (see "Moving the baseline").

The baseline is vg `202fc35d3`, the last commit that Benedict Paten regression-tested on
his own data. On both contigs, `caa925e81` and `f5b4e6b5c` reproduce it exactly. The baseline
`calls.vcf` comes from `21ce1a951`, which groups the INFO, FORMAT and FILTER lines of the VCF
header by type. Its records are the same as those of `202fc35d3`.

## Running it

Everything is a script in this directory, and a person can run them by hand exactly as
here. `config.sh` names every input and setting; override any of them from the
environment. First set `REGRESS_DIR` to a directory with about 5 GB free per contig for the gate's
data, tools and runs. It has no default. Then, from this directory:

1. **Install the scoring tools, once:** `./setup_tools.sh`. It installs truvari 5.5.0 and
   whatshap 2.8 in a venv at `$REGRESS_DIR/tools/venv` (`TOOLS_VENV`). It also needs
   `mafft`, which `truvari refine` runs, and GNU time at `/usr/bin/time`, and installs them
   with `apt-get` through `sudo` if they are missing. Put aardvark 1.0.0 on `PATH`, or set
   `AARDVARK` (see "The inputs" for where a copy was).
2. **Prepare the inputs, once:** `./prep.sh <any vg>`. It checks the input files' md5s,
   downloads and checks the truth, cuts each graph's component for the contig, extracts the
   reads, and slices the truth. It takes about 5 minutes for chr20 and 12 for chr6, and a
   rerun skips what is done. Run it once per contig.
3. **Build the vg to test:** `./build_vg.sh <commit>`, which prints the binary's path last.
   It clones the repository this skill is in, and builds in a new `mktemp -d` directory
   under `$TMPDIR` or `/tmp`. Keep `TMPDIR` on local disk, because a network filesystem is
   slow for a source tree. Delete the directory when done. The build takes about 30
   minutes at `JOBS=64`. To test uncommitted work, build the working tree as usual and
   pass its `bin/vg` instead.
4. **Run the gate:** `./gate.sh <vg> <label>`. Both arms run at once, at `THREADS=32`
   each, and need about 60 GiB of memory between them (70 GiB on chr6). On 64 free cores
   the gate takes about 6 minutes on chr20. On chr6 the calls take about 8 minutes and the
   scoring about 20, most of it `truvari refine`. Outputs go to `$REGRESS_DIR/runs/<label>/{short,ont}/`. It passes,
   exiting 0, when it prints "identical to the baseline" for all four items of each arm:
   `calls.vcf`, `mosaic.tsv`, `anchors.tsv` and the scores.

The gate's steps can also be run one at a time: `call.sh <vg> <arm> <dir>`, then
`score.sh <dir>`, then `compare.sh <dir> <arm>`. `score.sh` skips truvari or whatshap if
it cannot find them, and `compare.sh` then reports that the scores differ.

## The inputs

The `vgcall-wg-benchmark-wf` Snakemake workflow of the user `shnegi` made these. When the gate
was made, the workflow was at
`/private/groups/migalab/shnegi/vg_anchors_project/vgcall-wg-benchmark-wf` on the UCSC
Genomics Institute's cluster, and `config.sh` defaults to the files there. If they have
moved, find them by name or md5, and set `SHORT_GBZ`, `SHORT_GAF`, `ONT_GBZ` and `ONT_GAF`
(each GAF's `.tbi` must sit beside it). The graph is
Glenn Hickey's Oct 9 HPRC `pgrealign` gref graph, haplotype-sampled to 32 haplotypes by
each read set's own k-mers, with HG002 banned and CHM13 kept as the reference. Each read set
was mapped with Giraffe to its own sampled graph, then sorted and tabix-indexed. A GAF only
works with the GBZ it was mapped to. The GAFs have no `@RN` header, so their md5s are what
ties them to their graphs.

| file (under `.../vgcall-wg-benchmark-wf/results/`) | bytes | md5 |
|---|---|---|
| `sampled/novaseq-pcrfree-30x/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gbz` | 3,170,289,856 | `ddcd39170224632325fbb92d585b4f9e` |
| `aln/novaseq-pcrfree-30x/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gir-default.sorted.gaf.gz` | 24,287,745,993 | `db7028f3c1f434f75b5a236d00e7b717` |
| `.tbi` of that | | `c6705e9df8d03189cad0de28221e85f2` |
| `sampled/E821/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gbz` | 3,178,635,096 | `95fa448761c4025730a406ec8bd4f1c0` |
| `aln/E821/oct9-pgrealign-gref.hs32-refCHM13-banHG002.gir-r10.sorted.gaf.gz` | 67,403,338,564 | `2f3596942072ad41effbc45824c70dc7` |
| `.tbi` of that | | `a4ec37d5394ee5dbd741853215468eb5` |

`prep.sh` checks the GBZs and the indexes every time. Set `CHECK_GAFS=1` to check the GAFs
themselves too; that takes several minutes.

The baseline was scored with aardvark 1.0.0, the release binary
(`aardvark-v1.0.0-x86_64-unknown-linux-gnu`). `shnegi` kept a copy at
`/private/groups/migalab/shnegi/apps/aardvark/v1.0.0/aardvark-v1.0.0-x86_64-unknown-linux-gnu/aardvark`,
which `config.sh` falls back to if `aardvark` is not on `PATH`.

The truth is GIAB's HG002 v5.0q release on CHM13v2.0, `smvar` and `stvar`, each a VCF and a
benchmark BED. `prep.sh` downloads them and checks them against the md5s in the release's
`checksum.md5`.

### What prep.sh derives

For each arm and contig, under `$REGRESS_DIR/data/<arm>/`:

- `<contig>_0_<contig>.gbz`, from `vg chunk --gbz --contig <contig>`: the graph component
  that holds the contig.
- `<contig>.gaf.gz`: the GAF records whose every node is in that component.
- `<contig>.fa`: CHM13's sequence for the contig, from the graph. It must be the same in
  both graphs.

| contig | arm | component nodes (ID range) | GAF records | GAF content md5 (`zcat \| md5sum`) | FASTA md5 |
|---|---|---|---|---|---|
| chr20 | short | 2,551,512 (110548802–116688574) | 13,270,340 | `4388a9fbf48fccd7ccd759b9b6a5a25d` | `ac76aa812eba8f86cb8f0f80649dea4b` |
| chr20 | ont | 2,561,176 (110548802–116688574) | 85,605 | `f6aee1c98f87642d7e8ee58862f6013d` | `ac76aa812eba8f86cb8f0f80649dea4b` |
| chr6 | short | 5,507,411 (164206511–173206682) | 34,051,545 | `1de7b83a1c19a0625635655f7fb72076` | `0c3d3f0392ee84b55a83233f306e9938` |
| chr6 | ont | 5,514,432 (164206511–173206682) | 236,349 | `bc25caadf04925be1191ae772fc06142` | `0c3d3f0392ee84b55a83233f306e9938` |

The reads must be extracted. Pointing `vg call` at a contig's subgraph with the
whole-genome GAF makes it abort on an assertion in `gaf_to_alignment`
(`cur_offset < cur_len ...`), because some reads' paths leave the subgraph: on chr20, 163
short-read records and 2 ONT records. `extract_reads.py` drops those reads. On chr6 no
read leaves.

## The arms

Each arm is
`vg call <contig>_0_<contig>.gbz -p CHM13#0#<contig> -s HG002 -t 32 <arm options> --gaf-reads <contig>.gaf.gz --gaf-index <contig>.gaf.gz.tbi --anchors-out anchors.tsv --mosaic-out mosaic.tsv`:

- **short**: `--read-likelihood --mismap-max 0.95 --read-phasing --regenotype`, which is
  Benedict Paten's Illumina gate arm.
- **ont**: `--read-likelihood --preset ont`, as the `vgcall-wg-benchmark-wf` workflow runs it on
  the whole genome.

Both write anchors, so that the gate covers the anchor code on both read types.

## What the baseline must hit

The md5s are in `expected.tsv`. They were taken from `202fc35d3`, and `caa925e81` and
`f5b4e6b5c` reproduce them. The md5s of `calls.vcf` were taken from `21ce1a951`, whose header
has the same lines in another order. The headers of the mosaic and anchors files name the vg version and
the paths of the graph and reads, so their `#vg-version`, `#graph` and `#reads` lines are left out
of the md5s:

| contig | arm | `calls.vcf` | `mosaic.tsv` | `anchors.tsv` |
|---|---|---|---|---|
| chr20 | short | `42b225daaddc854d38641707fc200de9` (110,945 records) | `6ad75227e7dfa6a56adda5e37adf6d38` | `f27a6fcc5334e116dd3c72a0ac1101a2` |
| chr20 | ont | `ca4caf175f42f458b6ba4143b5678689` (118,127 records) | `8cb0405ab11fcb465caa3c860dc09af6` | `2f72b2671a1c823dc78732829c5e0de3` |
| chr6 | short | `d1af6962f9bbd4e293640894b96ebe49` (298,898 records) | `e4c05278623836422f252723e2c52e45` | `4ce059ee0f723c6063af1a37d524f698` |
| chr6 | ont | `b1cffcbcb42527029e9a9c1178fe5dfa` (300,803 records) | `51318e5d1e0e3a0f4c52e53087d7dd00` | `d83b5d69edb13f36bbe99a70b40b93c4` |

The output does not depend on the run or on the build. On chr20, a second run, a run at
`THREADS=16`, and a static build of `202fc35d3` by `shnegi` (built with `-march=nehalem`, and
kept at `/private/groups/migalab/shnegi/apps/vg-builds/vg_202fc35d3`) all gave the same
md5s as `202fc35d3`. The chr6 baseline was called at `THREADS=24`. So any difference comes from the
code. stderr matches too, once its lines are sorted, timings
are masked, and the output paths it names are ignored.

The scores are in `expected-score-<contig>-<arm>.json`. Inside the truth's benchmark
regions on each contig, they are:

| contig | arm | small F1 | SNV F1 | indel F1 | SV F1 | SV F1 after `truvari refine` | phasing: switches / flips (of assessed pairs) |
|---|---|---|---|---|---|---|---|
| chr20 | short | 0.979100 | 0.993138 | 0.928698 | 0.646108 | 0.720307 | 1,316 / 293 (61,359) |
| chr20 | ont | 0.974661 | 0.995645 | 0.898381 | 0.686825 | 0.832536 | 71 / 211 (63,480) |
| chr6 | short | 0.983179 | 0.994816 | 0.939860 | 0.715348 | 0.794446 | 2,623 / 543 (168,259) |
| chr6 | ont | 0.979159 | 0.996632 | 0.913476 | 0.743778 | 0.891423 | 11 / 336 (170,026) |

The small-variant scores are aardvark's genotype-level (`GT`) F1, for `ALL`, `Snv` and
`JointIndel`. Before scoring, records are split per ALT and kept only if PASS. Records under
50 bp are scored against `smvar` by aardvark, and all records against `stvar` by `truvari
bench --pick ac --passonly`, then by `truvari refine`. This is how the `vgcall-wg-benchmark-wf` workflow scores. The
phasing counts are `whatshap compare` against the truth VCF's phase over the whole contig. A flip
is one site phased against both of its neighbours, which whatshap's raw count reports as two
switch errors. The scoring is deterministic: identical calls gave identical scores.

As a rough guide, on chr20 one arm takes 3.5 to 4.5 minutes of wall time at 32 threads
with both arms running. The short arm peaks at about 27 GiB and the ONT arm at about 32 GiB
(15 and 19 GiB at 16 threads). On chr6 at 24 threads, with six calls sharing the machine,
the short arm took 6.5 to 7.5 minutes and peaked at about 35 GiB, and the ONT arm took
about 8 minutes and peaked at about 30 GiB. Timings taken while other jobs share the machine vary by 20% or more, so a
performance claim needs back-to-back runs with nothing else running. Run alone on a
160-core machine, the short arm took 3:10 to 3:37 of wall time and about 700 s of user
CPU, for both `202fc35d3` and `f5b4e6b5c`.

## Moving the baseline

When a change is meant to move the output, show that its scores are no worse, or that a
loss is a trade worth making, before it becomes the baseline. Then:

```
for arm in short ont; do cp $REGRESS_DIR/runs/<label>/$arm/score/summary.json expected-score-$CONTIG-$arm.json; done
```

Do this for each contig. Rewrite that contig's lines of `expected.tsv` from the run's files,
hashing `mosaic.tsv` and `anchors.tsv` without their `#vg-version`, `#graph` and `#reads` lines,
as `compare.sh` does. Then update the tables above, naming the new
baseline commit and what moved.

## Benedict Paten's data

Benedict's tier-2 gate (the private `benedictpaten/vg-call-eval` repository, `scripts/tier2/`)
used other files, which were not available when this gate was built. In PR 4990's
discussion (https://github.com/vgteam/vg/pull/4990#issuecomment-6301317554) he offered to
copy them wherever they were wanted, so look for them by the names and md5s below, or ask
him. If they turn up, the
same scripts can run on them. Set `SHORT_GBZ`, `SHORT_GAF`, `ONT_GBZ` and `ONT_GAF`, remove
the md5 checks in `prep.sh` (or replace them with the md5s below), and use a fresh
`REGRESS_DIR`. A GAF that is not sorted must first be sorted and indexed:

```
zcat X.gaf.gz | grep -v '^@' | awk -F'\t' '$6 != "*"' | vg gamsort -G - -t 8 | bgzip -@ 8 > X.sorted.gaf.gz
tabix -p gaf X.sorted.gaf.gz
```

The new outputs will not match `expected.tsv`, so first run `202fc35d3` on them to make a
baseline of their own.

His short-read tier-2 inputs, from PR 4990's discussion:

| file | bytes | md5 |
|---|---|---|
| `hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz` | 3,337,312,224 | `c44be569d4bebcc01f3b289da7451cc1` |
| `hprc-v2.1-mc-chm13-eval.HG002.hap32.gaf.gz` (not sorted; `@RN 987b3072…e923`) | 28,834,979,292 | `2cde3244c1f70fb874058535d50225e5` |
| `hprc-v2.1-mc-chm13-eval.HG002.gbz` (4 haplotypes) | 3,090,019,016 | `24b60d23b757ee77a5cdd964d875ff93` |
| `hprc-v2.1-mc-chm13-eval.HG002.gaf.gz` (`@RN 472bc25a…b145`) | 28,787,602,890 | `5b5ee214ea3afdccfe77401503efb3d9` |

His truth was the T2T-Q100 draft (`NIST_HG002_DraftBenchmark_defrabbV0.019-20241113`,
`CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz`, md5 `fdc6d410b904779e8f116b60e54f8e27`, and its
BED, md5 `6cfff916d4652a9c630700ea4424b965`). The v5.0q truth differs from it only slightly.

For `caa925e81` against `202fc35d3`, he reported these chr20 numbers. They come from "the
chr20 graph we used to gate `202fc35`", which is not the tier-2 34-haplotype graph. So they
are a reference, not a target:

| reads | ALL F1 | SNV F1 | INDEL F1 | SV F1 | SV refined F1 | switches / flips |
|---|---|---|---|---|---|---|
| ONT (`--preset ont`) | 0.979329 | 0.996218 | 0.917635 | 0.686118 | 0.831745 | 53 / 187 |
| Illumina (`--mismap-max 0.95 --read-phasing --regenotype`) | 0.981812 | 0.994688 | 0.936061 | 0.640189 | 0.705400 | 1,316 / 297 |

## What the gate does not cover

- `--gaf-base`. The gate reads through the tabix index. Commits `d9b5d874d` to `807474092`
  changed only the `--gaf-base` route, so this gate cannot see whether they changed output.
  Testing them needs `gaf-base` from jltsiren's gbz-base on `PATH`, and a database built with
  `gaf-base construct --preset long` (for ONT) from the extracted `<contig>.gaf.gz`.
- Contigs other than chr20 and chr6, other ploidies (`--ploidy-bed` for chrX and chrY), the support-based
  callers, and `-S gref_CHM13`. vg's own `make test` covers the support-based callers on
  small graphs.
