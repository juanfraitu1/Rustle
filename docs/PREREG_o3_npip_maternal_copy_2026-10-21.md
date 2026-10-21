# Pre-registration: does the reference-free chain name KB3781's maternal-only NPIP locus as exactly one new copy, blind to his assemblies? (2026-10-21)

**User (2026-10-21):** "ok lets do that". The proposal: cluster his NPIP reads, link alleles at 0.96%, and check that exactly one new
candidate appears and that it holds the 29 reads of register 1360. His assemblies serve only as the answer key.

**Question (O3: detect and flag).** Run on KB3781's own fibroblast reads of the NPIP family, with nothing masked, does the frozen
reference-free chain:
- name exactly one candidate new copy, namely the maternal-only locus whose reads minimap2 forces onto reference NPIP copies;
- hold that locus's reads;
- raise no false flag?

## The chain (frozen)

The chain is that of Amendments 7-9 of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`, with code in
`bench/rna_allele/link_test.py`, `merge_test.py` and `control_test.py`. It never sees his haplotype assemblies.

1. **Net.** Every read of his fibroblast library (`fibroblasts/GCA_029281585.2_flnc_mm.bam`, aligned to GGO with `-N 50 -p 0.1`)
   with any alignment record on one of the 25 NPIP copies (catalog family `NPIP`). Each read's sequence comes from its primary record,
   in sequencing orientation. Two deviations from Amendments 7-9:
   - **No 1,000-read cap.** The cap only bounded runtime: IsoCon ran the 2,951-read NPIP net in 93 s on 10-01, and a cap would
     subsample the target's reads.
   - **No unmapped reads.** His 959 unmapped reads have a median length of 69 bp, and 10-01 found no family evidence in them.
2. **IsoCon** 0.3.3: `IsoCon pipeline --nr_cores 4`, defaults.
3. **Flag.** An output is "not in the reference" when its best spliced alignment to GGO (`minimap2 -c -x splice:hq -uf -N 20`, the
   GGO splice index) has identity × coverage < 0.999.
4. **Link.** A flagged output whose whole-length divergence to that best locus is d <= δ = 0.00958 is an allele of it; otherwise it
   is a new-copy contig. δ is the 99th percentile of exonic divergence between his two haplotypes over 28,541 single-copy genes,
   frozen on 10-01.
5. **Merge.** New-copy contigs are aligned all-vs-all (`minimap2 -c -x asm20 --cs -N 200 -p 0.1`). Two join when their best
   alignment covers >= 50% of the shorter contig and de <= δ. The components are the candidate copies.
6. **Support floor.** A candidate counts only with >= 2 contigs. Amendment 9's no-deletion control proposed this floor, after the
   fact, for the next tests; it is fixed here before any run. One-contig candidates are listed, not counted.
7. **A candidate's reads:** the reads IsoCon used for its contigs (`cluster_info.tsv`).

## Answer key (his assemblies; evaluation only)

- **Candidate class** (Amendment 9's `classify`): the best contig's spliced alignment to his maternal and paternal assemblies
  (`-c -x splice:hq -uf -N 20`, the haplotype splice indexes).
  - identity × coverage < 0.999: **unmatched** (false);
  - on the haplotype the reference is identical to: **pri** (false);
  - inside the lift of an NPIP copy (the frozen 10-01 asm5 lift, `truth_lift.py`, copy spans): **allele** (false);
  - otherwise: **haplotype-only** (true).
- **Target T:** the maternal-only locus of register 1360, his maternal chr18 `CM054600.2:6,720,000-6,820,000`.
- **Target reads R_T:** net reads whose best alignment in his genome lies in T and scores >= 5 `ms` above their best alignment
  outside T. Alignments are `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1` to each haplotype assembly, compared by `ms`;
  alignments with an intron of >= 100 kb are dropped. The 29 tied reads of register 1360 are reported among them.

## Bars

- **E1 (found):** at least one counted candidate is haplotype-only, with its hit inside T.
- **E2 (one copy, one candidate):** exactly one such candidate.
- **E3 (holds its reads):** that candidate holds at least half of R_T.
- **E4 (pure):** at least 90% of its reads are in R_T.
- **E5 (no false flag):** no other counted candidate is false (unmatched, pri or allele). Other haplotype-only candidates are
  reported as additional flags, not failures; one example would be 10-01's maternal-only copy at chr18 mat 17.40 Mb.

**PASS** = E1-E5. **FAIL** otherwise, and the report names the bars that fail.

## Seen before: this is a positive control of the chain, not a blind discovery

- Register 1360 found the locus from his assemblies: 29 tied NPIP reads with 0-3 substitutions there and 18-167 at their best
  reference placement.
- On 10-01, IsoCon on the same library (with NPIPA2 masked) produced 4 outputs from this sequence (chr18 mat 6.72-6.76 Mb, 79-93% to
  the reference), and Amendment 5 moved 31 reads onto two of them (iso_31 / iso_32).
- **Never run before:**
  - the full chain (flag, link, merge, support floor) on NPIP with nothing masked;
  - its read accounting against R_T;
  - its false flags within NPIP.

## Rules

- Run once; bugs → amendment + rerun.
- Heavy commands in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Gorilla only. Push only when asked.
- Code: a new driver, `bench/rna_allele/npip_o3_test.py`, calls the frozen pieces (the flag and link of `control_test.contigs`,
  `merge_test.pairs` / `components`, `truth_lift`), with tests for its own new parts (the net, R_T and the bars).
