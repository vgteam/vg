# Machine-written write-ups

These write-ups of the read-likelihood caller were written by a coding agent and have not been
read and cleared by a person. See `artifacts/README.md`.

## Method write-ups

Long-form descriptions of methods live here in the source tree, beside the code that implements
them, so that a method and its description are reviewed in one diff and versioned together. A
write-up covers what a user needs to understand and configure a method: the model it uses, how
each part is computed, the options that control each part, and what its outputs mean. Default
values live in the code and the command's `--help`, not in the write-up. Rationale that only
matters to someone editing the code stays in the headers.

- [read-likelihood-genotyping.md](read-likelihood-genotyping.md) — the `vg call --read-likelihood`
  caller as a whole, and the place to start: its vocabulary, how its passes and rounds are
  ordered, how a site's genotype is settled, nested sites, phasing, the output fields and files,
  and the options for each part.
- [read-likelihood-direct-genotyping.md](read-likelihood-direct-genotyping.md) — how
  `--read-likelihood` genotypes a site from its reads: the site likelihood and how each of its
  terms is computed.
- [read-likelihood-linkage-model.md](read-likelihood-linkage-model.md) — the linkage model, which
  re-decides genotypes and phases them from the panel haplotypes, and how it differs from
  PanGenie's.
- [read-likelihood-read-phasing.md](read-likelihood-read-phasing.md) — read phasing, which
  re-decides the phases from reads that span several sites, and re-genotyping, which uses the phase
  to correct each site's likelihoods.

The [wiki](https://github.com/vgteam/vg/wiki) remains the home for tutorials and worked examples.
It is a separate repository (mounted here as the `wiki` submodule), so anything published there
cannot be reviewed alongside a code change; a page there should link to the write-up rather than
restate it.

## Code structure

Developer-facing, and the exception to the rule above: a description of how a subsystem's files
depend on each other is too long for a header and belongs to no single one of them.

- [read-likelihood-architecture.md](read-likelihood-architecture.md) — what the
  `--read-likelihood` files are, the dependency graph between them, and where the coupling is.
