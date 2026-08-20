# Diploid-aware surjection: implementation and evaluation summary

## Summary

`vg surject` converts an alignment against a graph into an alignment against one or more linear paths in that graph. Ordinary surjection evaluates each input graph alignment against the requested paths and reports the best linear projection. This becomes difficult when the graph contains two haplotypes for the target sample and Giraffe has retained several possible graph mappings for a read.

`vg surject --diploid-map SAMPLE` addresses this by considering the two sample haplotypes together. It groups all graph mappings belonging to a read, projects them onto the sample's haplotype paths, compares the resulting alternatives globally, and chooses a single overall primary alignment. For paired short reads, it also evaluates the two mates jointly and can learn the fragment-length distribution from the data.

Two complementary evaluations were performed:

1. **HG002 acceptance testing:** real HiFi and Illumina reads were mapped, surjected to CHM13, called with DeepVariant, and compared with GIAB truth variants. These results measure downstream variant-calling effects.
2. **KOLF2.1J simulation:** reads were simulated from a known linear KOLF2.1J assembly, mapped to a graph containing KOLF2.1J haplotypes and CHM13, and compared directly with their known source positions after surjection. These results measure read-placement accuracy and MAPQ calibration.

The non-diploid changes were nearly neutral in direct KOLF placement accuracy. Diploid-map produced a small net gain for HiFi and a much larger gain for paired short reads, where joint mate selection and fragment-length information are especially useful.

## What `--diploid-map` does

The command-line form is:

```bash
vg surject \
    -x graph.gbz \
    --diploid-map SAMPLE \
    --bam-output \
    input.gam
```

Conceptually, `--diploid-map SAMPLE` behaves like selecting the sample assembly with `--into-ref SAMPLE`, but it additionally enables multimap-aware, haplotype-aware selection.

For each read it:

1. Collects the primary and alternative graph mappings sharing that read name.
2. Surjects every graph mapping onto the relevant paths from both sample haplotypes.
3. Compares the haplotype projections for each graph mapping.
4. Compares all projected alternatives globally.
5. Emits one overall primary alignment and marks the remaining overlapping alternatives secondary.
6. Computes mapping confidence from competition among the projected alternatives.
7. Records haplotype and alignment-quality information in `hp`, `hq`, and `aq` tags.

The tags distinguish three related quantities:

- `hp`: whether a projection is the preferred or alternative haplotype projection.
- `hq`: confidence in the haplotype choice.
- `aq`: the mapping quality of the original graph alignment before diploid-aware projection.

### Paired short reads

With interleaved paired input, diploid-map evaluates candidate placements for the two mates jointly. It learns a fragment-length distribution from high-confidence pairs and adds a fragment-length likelihood to the pair score. This helps distinguish haplotype/path combinations that have similar single-read alignment scores but differ in whether they form a plausible fragment.

The paired implementation buffers ambiguous fragments while the distribution is being learned, then replays them after an estimate is available. If a reliable distribution cannot be learned, it falls back without the fragment-length contribution.

## What the non-diploid changes are

The candidate vg image contains several `surject` improvements that apply even when `--diploid-map` is not enabled. They should be separated conceptually from diploid-aware selection.

### Mapper-declared tail pruning

The new `--prune-tail-region` option removes surjection anchors lying entirely inside tail regions declared by the mapper. Long-read tails can contain weak or repetitive anchors that pull the projected alignment toward an incorrect local placement. Removing these anchors allows the tail to be realigned against the selected linear path.

### Better alternative-alignment classification

When several path projections exist, alternatives overlapping the same query interval are now classified as secondary. Disjoint query-interval projections can be retained as supplementary alignments. This is more faithful than treating all non-primary projections as the same kind of alternative.

### Path ranking

Candidate paths are ranked using their best available projected alignment rather than an ordering that can be influenced by weaker alternatives. This makes primary selection less sensitive to path enumeration order.

### Substitution-anchor realignment

Anchors containing substitutions are realigned instead of being accepted as exact fixed anchors. This reduces cases in which an imperfect anchor constrains the rest of the surjection incorrectly.

### Target-path-aware anchor sliding

The anchor-pruning logic now checks whether an anchor sequence occurs again at a nearby offset on the target path, not only whether a similar shift is possible in the read sequence. This more directly detects ambiguous repeat anchors on the path being used for BAM output.

Adam's acceptance-test candidate tested these non-diploid changes with ordinary surjection and `--prune-tail-region`; it did not enable `--diploid-map`.

## HG002 acceptance tests: Adam's non-diploid candidate

The acceptance workflow uses real HG002 reads. It maps them to the HPRC graph, surjects to CHM13, calls variants with DeepVariant, and evaluates the VCF against GIAB truth with Aardvark. Therefore, FN and FP below refer to **variant calls**, not individual read placements.

### HiFi

| Metric | WDL baseline | Non-diploid candidate | Change |
|---|---:|---:|---:|
| False negatives | 47,607 | 47,702 | +95 |
| False positives | 29,578 | 29,214 | -364 |
| Recall | 98.8851% | 98.8829% | -0.0022 pp |
| Precision | 99.3049% | 99.3134% | +0.0085 pp |
| F1 | 99.0946% | 99.0977% | +0.0031 pp |
| Genotype FN | 7,596 | 7,599 | +3 |
| Genotype FP | 3,234 | 3,269 | +35 |

The HiFi candidate trades 95 additional false negatives for 364 fewer false positives. The net F1 change is positive but very small.

### Illumina

| Metric | WDL baseline | Non-diploid candidate | Change |
|---|---:|---:|---:|
| False negatives | 74,653 | 74,604 | -49 |
| False positives | 15,818 | 15,748 | -70 |
| Recall | 98.2518% | 98.2529% | +0.0011 pp |
| Precision | 99.6245% | 99.6262% | +0.0017 pp |
| F1 | 98.9334% | 98.9348% | +0.0014 pp |
| Genotype FN | 6,065 | 6,038 | -27 |
| Genotype FP | 5,344 | 5,362 | +18 |

The non-diploid Illumina candidate slightly improves FN, FP, recall, precision, and F1. As with HiFi, the changes are small, which is expected for a targeted surjection change evaluated after a full variant-calling pipeline.

## HG002 short-read acceptance test with `--diploid-map`

The short-read diploid candidate uses the same overall acceptance workflow but enables diploid-aware CHM13 surjection. The paired implementation includes fragment-distribution learning.

| Metric | Baseline in diploid run | Diploid-map candidate | Change |
|---|---:|---:|---:|
| False negatives | 74,648 | 74,602 | -46 |
| False positives | 15,818 | 15,745 | -73 |
| Recall | 98.2519% | 98.2530% | +0.0011 pp |
| Precision | 99.6245% | 99.6263% | +0.0017 pp |
| F1 | 98.9335% | 98.9348% | +0.0014 pp |
| Genotype FN | 6,065 | 6,039 | -26 |
| Genotype FP | 5,343 | 5,360 | +17 |

The diploid-map result is almost identical to Adam's non-diploid candidate at the VCF level: it has two fewer FN, three fewer FP, and an F1 difference of only about 0.00006 percentage points. Thus, the acceptance test shows no regression and a tiny numerical improvement, but it does not by itself expose the much larger read-selection changes seen in controlled simulation.

The completed short-read acceptance comparison did not use the later problematic HiFi MM5 baseline configuration. Its two baseline results agree within five FN and one genotype FP.

## KOLF2.1J calibration experiment

### Purpose

The KOLF experiment isolates surjection more directly than the acceptance test. A diploid graph contains the two KOLF2.1J haplotypes and CHM13. Reads are simulated from a separate linear KOLF2.1J assembly, giving every read a known source position. The simulated reads are mapped back to the graph with Giraffe while retaining multiple graph mappings. Surjected primary alignments are then compared with truth using:

```bash
vg gamcompare -T -r 50 surjected.gam truth.gam
```

A read is considered correct when its surjected position is within 50 bp of its truth position.

Three modes are compared:

1. **WDL baseline:** ordinary surjection with `--into-ref KOLF2.1J --multimap`.
2. **Non-diploid candidate:** candidate code with the same non-diploid target and multimap behavior.
3. **Diploid candidate:** candidate code with `--diploid-map KOLF2.1J`.

The clean WDL baseline image identifies itself as `vg v1.76.1-114-g431e126`. Earlier directories were named `v175`, but that label is inaccurate; the exact commit/version file should be used when describing the control.

## KOLF HiFi comparison

One million simulated HiFi reads were evaluated.

| Condition | Correct | Incorrect | Accuracy | Net versus baseline |
|---|---:|---:|---:|---:|
| WDL baseline | 966,204 | 33,796 | 96.6204% | — |
| Non-diploid candidate | 966,194 | 33,806 | 96.6194% | -10 |
| Diploid-map candidate | **966,628** | **33,372** | **96.6628%** | **+424** |

### Rescued and regressed reads

| Comparison | Rescued | Regressed | Net |
|---|---:|---:|---:|
| Baseline to non-diploid | 11 | 21 | -10 |
| Non-diploid to diploid-map | 16,078 | 15,644 | +434 |
| Baseline to diploid-map | 16,078 | 15,654 | +424 |

For HiFi, the non-diploid changes are essentially neutral for truth-position accuracy. Diploid-map changes the selected placement for more than 31,000 reads, but rescues only 424 more reads than it regresses. The overall accuracy gain is therefore small: 0.0424 percentage points.

This does not mean diploid-map has little effect. The large rescue and regression counts show that it frequently selects a different haplotype/path placement; the favorable and unfavorable changes nearly balance for these long reads.

## KOLF paired short-read comparison

The short-read evaluation contains one million reads, or 500,000 paired fragments.

| Condition | Correct | Incorrect | Accuracy |
|---|---:|---:|---:|
| Clean WDL baseline | 556,642 | 443,358 | 55.6642% |
| Historical non-diploid table | 556,804 / 999,599 rows | 442,795 | approximately 55.70% |
| Diploid-map candidate | **603,963** | **396,037** | **60.3963%** |

The historical non-diploid table has a provenance limitation: it contains 999,514 unique read IDs, 81 duplicated IDs contributing 85 extra rows, and is missing 486 baseline read IDs. For a fair three-way transition analysis, duplicated IDs are excluded and only IDs present exactly once in all three conditions are used.

### Matched three-way transitions

| Comparison | Rescued | Regressed | Net |
|---|---:|---:|---:|
| Baseline to non-diploid | 245 | 90 | +155 |
| Non-diploid to diploid-map | 226,440 | 179,444 | +46,996 |
| Baseline to diploid-map | 226,540 | 179,389 | +47,151 |

On the complete baseline and diploid tables, without involving the incomplete non-diploid table, diploid-map rescues 226,747 reads and regresses 179,426. The net gain is 47,321 correct reads, corresponding to:

- +4.7321 percentage points in read-placement accuracy.
- 10.67% fewer incorrectly placed reads.
- 51.14% of baseline errors rescued.
- 32.23% of baseline-correct reads regressed.

The bidirectional changes remain large, but unlike HiFi they do not nearly cancel. Paired diploid-map produces a substantial net benefit.

### Fragment-level result

| Fragment outcome | Baseline | Diploid-map | Change |
|---|---:|---:|---:|
| Both mates correct | 256,869 | **301,966** | **+45,097** |
| Exactly one mate correct | 42,904 | **31** | **-42,873** |
| Neither mate correct | 200,227 | **198,003** | **-2,224** |

This is the clearest evidence for the paired diploid method. Ordinary surjection often produces fragments with one correct mate and one incorrect mate. Diploid-map almost eliminates that mixed category by scoring the mates together. Most of those fragments become correct for both mates, although some become incorrect for both.

### Mapping-quality calibration

The baseline assigns MAPQ 60 to 903,660 reads, but only 57.68% of those alignments satisfy the 50 bp truth criterion. Diploid-map assigns MAPQ 60 to only 86,808 reads, of which 86,802 are correct (99.993%). It places most ambiguous reads in the MAPQ 1–9 range instead of reporting them as highly confident.

Mean MAPQ falls from 54.85 for baseline to 11.88 for diploid-map, while positional accuracy increases. The lower mean is therefore not a loss of quality; it reflects more conservative and better-separated confidence estimates.

## KOLF HiFi speed and memory benchmark

A separate benchmark measured surjection runtime and peak resident memory on the one-million-read KOLF HiFi mappings. It used the exact WDL baseline and candidate container images and compared:

1. `baseline_wdl`: WDL baseline image with `--into-ref KOLF2.1J --multimap`.
2. `candidate_control`: candidate image with the same non-diploid mode.
3. `candidate_diploid`: candidate image with `--diploid-map KOLF2.1J`.

Both GAM and GAF input were tested. Each condition was run five times. Condition order was alternated between replicates, and all conditions for a given input format ran on the same compute node. Every replicate completed successfully.

### GAF input

| Condition | Replicates | Mean elapsed time | Change versus baseline | Mean peak RSS |
|---|---:|---:|---:|---:|
| WDL baseline | 5 | 804.1 s (13.40 min) | — | 64.4 GiB |
| Non-diploid candidate | 5 | **724.6 s (12.08 min)** | **9.9% faster** | 65.4 GiB |
| Diploid-map candidate | 5 | **739.3 s (12.32 min)** | **8.1% faster** | 66.5 GiB |

For GAF, diploid-map is about 2.0% slower than the non-diploid candidate but remains 8.1% faster than the WDL baseline. Peak-memory differences are modest: diploid-map uses about 2.1 GiB more than baseline.

### GAM input

| Condition | Replicates | Mean elapsed time | Change versus baseline | Mean peak RSS |
|---|---:|---:|---:|---:|
| WDL baseline | 5 | 674.3 s (11.24 min) | — | 70.9 GiB |
| Non-diploid candidate | 5 | **595.6 s (9.93 min)** | **11.7% faster** | 69.3 GiB |
| Diploid-map candidate | 5 | **626.2 s (10.44 min)** | **7.1% faster** | **93.1 GiB** |

For GAM, diploid-map is about 5.1% slower than the non-diploid candidate but remains 7.1% faster than baseline. The important cost is memory: diploid-map uses about 22.2 GiB more than baseline and 23.8 GiB more than the non-diploid candidate. This is consistent with grouping multiple mappings per read and retaining candidate placements for global haplotype selection.

### Performance interpretation

The speed study does not show a runtime regression relative to the WDL baseline. The newer non-diploid code is fastest, while diploid-aware selection adds a modest overhead relative to that candidate control. The main resource concern is GAM-input peak memory. GAF input has a much smaller memory penalty and may be preferable when memory is constrained, subject to equivalent correctness and metadata handling.

The source benchmark tables are:

- `long/WDL_VG_SURJECT_BENCHMARK/gaf/summary.36802132.tsv`
- `long/WDL_VG_SURJECT_BENCHMARK/gam/summary.36802133.tsv`

## Overall interpretation

The experiments support four conclusions:

1. **The non-diploid changes are safe but have small aggregate effects.** Adam's acceptance tests show tiny F1 improvements, and the KOLF simulations show nearly neutral direct placement accuracy.
2. **Diploid-map changes primary selection much more substantially than the non-diploid fixes.** This is visible in the large rescued and regressed read sets.
3. **The strongest benefit is for paired short reads.** Joint haplotype and fragment scoring increases KOLF positional accuracy by 4.73 percentage points and almost eliminates fragments with only one correct mate.
4. **Confidence calibration improves markedly in the short-read simulation.** High MAPQ becomes restricted to a smaller set of almost uniformly correct placements.
5. **Runtime remains competitive with baseline.** Diploid-map is 7–8% faster than the WDL baseline in the controlled KOLF HiFi benchmark, although GAM-input peak memory increases to about 93 GiB.

The HG002 acceptance results remain important because they show that these large read-level changes do not cause a downstream variant-calling regression. The KOLF simulation explains what the acceptance metrics cannot: diploid-map is changing haplotype/path selection and pair consistency even when the final VCF difference is small.

## Limitations

- The HG002 acceptance test evaluates variants, not direct read-placement truth.
- KOLF correctness is defined as position within 50 bp; it does not independently assess base-level alignment or CIGAR accuracy.
- The historical short-read non-diploid output should be regenerated because of missing and duplicated read IDs.
- Exact candidate version files were not archived beside every historical KOLF output.
- The clean WDL baseline is `v1.76.1-114-g431e126`, not an exact vg 1.75 release.
- Simulated KOLF reads and a graph containing the correct KOLF haplotypes represent a favorable setting; real samples and graphs with incomplete haplotype representation may behave differently.
