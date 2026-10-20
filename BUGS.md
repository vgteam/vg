# Known bugs on this branch

Bugs found in the read-likelihood caller (PR 4990) and its merge of PR 5057, still to fix. Each
fix crosses its row off in the same commit. Delete this file once every row is fixed, before the
branch merges. Rows marked TODO have a matching `TODO:` comment in the code.

| # | Where | Problem | Source | Status |
|---|---|---|---|---|
| 1 | `src/site_read_source.cpp` `thread_state()` | Resizing `threads` can move a `ThreadState` that another thread is using. Likely a data race. | reads agent (TODO) | open |
| 2 | `src/site_read_source.cpp` `deliver()` | Hands a mapping over twice if the caller's ranges overlap. | reads agent (TODO) | open |
| 3 | `src/site_read_source.cpp` `for_each_read_start` | Answers from holds only when off-panel detection is on. Should it always? (Question, maybe not a bug.) | reads agent (TODO) | open |
| 4 | `src/linkage_model.cpp` `pair_viterbi` | Backpointers are `uint16_t`, so a panel of more than 65536 haplotypes wraps them, and nothing checks the panel size. | phasing/linkage agent (TODO) | open |
| 5 | `src/round_history.*` | Rounds are compared by a hash of the state, so a collision reports a cycle that is not one. Documented, not changed. | phasing agent | open |
| 6 | `vg call --read-likelihood -r <snarls>` | `SnarlManagerDecomposition` now throws on unary snarls, which the old adapter showed. A CactusSnarlFinder snarl file would now fail. | merge of PR 5057 | open |
| ~~7~~ | ~~PR 4990 own files: `mosaic_writer.cpp`, `linkage_model.cpp`, `read_likelihood_caller.cpp`, `off_panel.cpp`, `site_values.cpp`~~ | ~~Use `<format>`/`std::format` or `using enum`, which GCC 10 (vg's minimum) lacks. CI would fail as PR 5057 did.~~ | ~~grep after PR 5057's CI fix~~ | fixed |
| ~~8~~ | ~~merge branch~~ | ~~Merge `snarl-manager-decomposition` again to pick up its GCC 10 fix (daf7987e2).~~ | ~~PR 5057 CI~~ | fixed |
| 9 | `src/vcf_record.cpp` `get_ref_interval` | The last step of the path is never compared, so a site whose end bound is the last step is not found. | records agent (TODO) | open |
| 10 | `src/vcf_record.cpp` `get_ref_interval` | After one backward scan, later forward scans still look for the flipped end. | records agent (TODO) | open |
| 11 | `flatten_common_allele_ends` | A `len_override` longer than the shortest allele reads past its end. | records agent (TODO) | open |
| 12 | `src/block_records.cpp` `SiteAlignment` | A genotype entry past the end of `called_traversals` is skipped, but later code still reads that traversal. | records agent (TODO) | open |
| 13 | `src/mosaic_writer.cpp` `PanelWalker::position_key` | Keeps only 20 bits of the haplotype index, and `position_at` returns a cache hit without checking the haplotype, so with more than 2^20 haplotypes one haplotype gets another's position (crosstalk, not just contention). Latent for real panel sizes. Fix: skip the cache for haplotype indices that do not fit, as `sequence_at` does for offsets. | records agent (TODO) | open |
| 14 | `RecordRenderer::reader` | Never read; removing it and `configure` needs a change in `multipass_caller.cpp`. Clean-up, not a bug. | records agent (TODO) | open |
| 15 | `src/site_scheduler.hpp` `SiteCallTally` | The "nested snarls" progress line prints the top-level count. | nested agent (TODO) | open |
| 16 | `StagedSiteTable::hand_off` | `panel_cache.clear()` frees no memory, and leaves `panel_cached` true. | nested agent (TODO) | open |
| 17 | `chain_key_of` | A hash, so two chains can share a key. | nested agent (TODO) | open |
| 18 | `emit_variant` (vcf_output_caller) | Takes the reference offset as an `int`, while FlowCaller and NestedFlowCaller pass a `size_t`, so an offset of 2^31 or more comes out wrong. | callers agent (TODO) | open |
