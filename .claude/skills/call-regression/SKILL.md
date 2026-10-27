---
name: call-regression
description: Regression-test vg call's read-likelihood caller (PR 4990's --read-likelihood, nesting, phasing, anchors and mosaic) on real HG002 reads, short and ONT, on chr20, by checking that the VCF, mosaic and anchors are byte-identical to a baseline and by scoring them against GIAB v5.0q. Use before pushing or merging any change that could touch vg call's output, and to set a new baseline after an intended output change.
---

# Regression-test vg call on HG002 chr20

The gate calls HG002 on chr20 twice with the vg under test: once from NovaSeq 30x short
reads and once from E821 ONT reads. Each read set has its own sampled graph. The gate then
checks the output files against the baseline's md5s and scores the calls against GIAB's
HG002 v5.0q benchmark.

For a change that should not move the output, such as a refactor, a speed-up or a change to
another subcommand, every file must be byte-identical, and that is the whole test. For a
change that is meant to move the output, the scores say whether accuracy rose or fell. Once
the change is accepted, its run becomes the new baseline (see "Moving the baseline").

The baseline is vg `202fc35d3`, the last commit that Benedict Paten regression-tested on
his own data. On these inputs, `caa925e81` and `f5b4e6b5c` reproduce it exactly.

## Running it

Everything is a script in this directory, and a person can run them by hand exactly as
here. `config.sh` names every input and setting; override any of them from the
environment. First set `REGRESS_DIR` to a directory with about 5 GB free for the gate's
data, tools and runs. It has no default. Then, from this directory:

1. **Install the scoring tools, once:** `./setup_tools.sh`. It installs truvari 5.5.0 and
   whatshap 2.8 in a venv at `$REGRESS_DIR/tools/venv` (`TOOLS_VENV`). It also needs
   `mafft`, which `truvari refine` runs, and GNU time at `/usr/bin/time`, and installs them
   with `apt-get` through `sudo` if they are missing. Put aardvark 1.0.0 on `PATH`, or set
   `AARDVARK` (see "The inputs" for where a copy was).
2. **Prepare the inputs, once:** `./prep.sh <any vg>`. It checks the input files' md5s,
   downloads and checks the truth, cuts each graph's chr20 component, extracts the reads,
   and slices the truth. It takes about 5 minutes, and a rerun skips what is done.
3. **Build the vg to test:** `./build_vg.sh <commit>`, which prints the binary's path last.
   It clones the repository this skill is in, and builds in a new `mktemp -d` directory
   under `$TMPDIR` or `/tmp`. Keep `TMPDIR` on local disk, because a network filesystem is
   slow for a source tree. Delete the directory when done. The build takes about 30
   minutes at `JOBS=64`. To test uncommitted work, build the working tree as usual and
   pass its `bin/vg` instead.
4. **Run the gate:** `./gate.sh <vg> <label>`. Both arms run at once, at `THREADS=32`
   each, and need about 60 GiB of memory between them. The gate takes about 6 minutes on
   64 free cores. Outputs go to `$REGRESS_DIR/runs/<label>/{short,ont}/`. It passes,
   exiting 0, when it prints "identical to the baseline" for all four items of each arm:
   `calls.vcf`, `mosaic.tsv`, `anchors.tsv` and the scores.

The gate's steps can also be run one at a time: `call.sh <vg> <arm> <dir>`, then
`score.sh <dir>`, then `compare.sh <dir> <arm>`. `score.sh` skips truvari or whatshap if
it cannot find them, and `compare.sh` then reports that the scores differ.

## The inputs

Shlok Negi made these with his `vgcall-wg-benchmark-wf` Snakemake workflow. When the gate
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
(`aardvark-v1.0.0-x86_64-unknown-linux-gnu`). Shlok kept a copy at
`/private/groups/migalab/shnegi/apps/aardvark/v1.0.0/aardvark-v1.0.0-x86_64-unknown-linux-gnu/aardvark`,
which `config.sh` falls back to if `aardvark` is not on `PATH`.

The truth is GIAB's HG002 v5.0q release on CHM13v2.0, `smvar` and `stvar`, each a VCF and a
benchmark BED. `prep.sh` downloads them and checks them against the md5s in the release's
`checksum.md5`.

### What prep.sh derives

| file (under `$REGRESS_DIR/data/<arm>/`) | what it is | check |
|---|---|---|
| `chr20_0_chr20.gbz` | `vg chunk --gbz --contig chr20`: the graph component that holds chr20 | the short arm's has 2,551,512 nodes and the ONT arm's 2,561,176, both with IDs 110548802 to 116688574 |
| `chr20.gaf.gz` | the GAF records whose every node is in that component | content md5 (`zcat \| md5sum`): short `4388a9fbf48fccd7ccd759b9b6a5a25d`, 13,270,340 records; ONT `f6aee1c98f87642d7e8ee58862f6013d`, 85,605 records |
| `chr20.fa` | CHM13's chr20, from the graph | md5 `ac76aa812eba8f86cb8f0f80649dea4b`, the same in both graphs |

The reads must be extracted. Pointing `vg call` at a contig's subgraph with the
whole-genome GAF makes it abort on an assertion in `gaf_to_alignment`
(`cur_offset < cur_len ...`), because some reads' paths leave the subgraph: 163 short-read
records and 2 ONT records here. `extract_reads.py` drops those reads.

## The arms

Each arm is
`vg call chr20_0_chr20.gbz -p CHM13#0#chr20 -s HG002 -t 32 <arm options> --gaf-reads chr20.gaf.gz --gaf-index chr20.gaf.gz.tbi --anchors-out anchors.tsv --mosaic-out mosaic.tsv`:

- **short**: `--read-likelihood --mismap-max 0.95 --read-phasing --regenotype`, which is
  Benedict Paten's Illumina gate arm.
- **ont**: `--read-likelihood --preset ont`, as Shlok's whole-genome workflow runs it.

Both write anchors, so that the gate covers the anchor code on both read types.

## What the baseline must hit

The md5s are in `expected.tsv`. They were taken from `202fc35d3`, and `caa925e81` and
`f5b4e6b5c` reproduce them:

| arm | `calls.vcf` | `mosaic.tsv` | `anchors.tsv` without its `#vg-version` line |
|---|---|---|---|
| short | `28d7ce959bb875d5e6830fe2b3b09e41` (110,945 records) | `921d3dc80f660fa28333c1f7028317e9` | `3a9a868c13c02ea9e2e1345bffcfef22` |
| ont | `1e9dce357c2c42192685e80f17c03497` (118,127 records) | `f72d665575b71f6815d1614fe282e8ef` | `be2ed774a8c69f1889ff3b0c0502c53a` |

The output does not depend on the run or on the build. A second run, a run at
`THREADS=16`, and Shlok's static build of `202fc35d3` (built with `-march=nehalem`, and
kept at `/private/groups/migalab/shnegi/apps/vg-builds/vg_202fc35d3`) all gave these same
md5s. So
any difference comes from the code. stderr matches too, once its lines are sorted, timings
are masked, and the output paths it names are ignored.

The scores are in `expected-score-short.json` and `expected-score-ont.json`. Inside the
truth's benchmark regions on chr20, they are:

| arm | small F1 | SNV F1 | indel F1 | SV F1 | SV F1 after `truvari refine` | phasing: switches / flips (of assessed pairs) |
|---|---|---|---|---|---|---|
| short | 0.979100 | 0.993138 | 0.928698 | 0.646108 | 0.720307 | 1,316 / 293 (61,359) |
| ont | 0.974661 | 0.995645 | 0.898381 | 0.686825 | 0.832536 | 71 / 211 (63,480) |

The small-variant scores are aardvark's genotype-level (`GT`) F1, for `ALL`, `Snv` and
`JointIndel`. Before scoring, records are split per ALT and kept only if PASS. Records under
50 bp are scored against `smvar` by aardvark, and all records against `stvar` by `truvari
bench --pick ac --passonly`, then by `truvari refine`. This is how Shlok's workflow scores. The
phasing counts are `whatshap compare` against the truth VCF's phase over all of chr20. A flip
is one site phased against both of its neighbours, which whatshap's raw count reports as two
switch errors. The scoring is deterministic: identical calls gave identical scores.

As a rough guide, one arm takes 3.5 to 4.5 minutes of wall time at 32 threads with both arms
running. The short arm peaks at about 27 GiB and the ONT arm at about 32 GiB (15 and 19 GiB
at 16 threads). Timings taken while other jobs share the machine vary by 20% or more, so a
performance claim needs back-to-back runs with nothing else running. Run alone on a
160-core machine, the short arm took 3:10 to 3:37 of wall time and about 700 s of user
CPU, for both `202fc35d3` and `f5b4e6b5c`.

## Moving the baseline

When a change is meant to move the output, show that its scores are no worse, or that a
loss is a trade worth making, before it becomes the baseline. Then:

```
for arm in short ont; do cp $REGRESS_DIR/runs/<label>/$arm/score/summary.json expected-score-$arm.json; done
```

Rewrite `expected.tsv` from that run's files, hashing `anchors.tsv` without its
`#vg-version` line as `compare.sh` does. Then update the tables above, naming the new
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
BED, md5 `6cfff916d4652a9c630700ea4424b965`). Shlok's v5.0q differs from it only slightly.

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
  `gaf-base construct --preset long` (for ONT) from the extracted `chr20.gaf.gz`.
- Other contigs, other ploidies (`--ploidy-bed` for chrX and chrY), the support-based
  callers, and `-S gref_CHM13`. vg's own `make test` covers the support-based callers on
  small graphs.
