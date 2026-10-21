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

## Results (driver a7a81ece, `bench/rna_allele/npip_o3.py`; run once; `/mnt/linuxdisk/tmp/o3_npip/`): **FAIL** (E3, E4, E5)

- **Chain:**
  - net 2,955 reads (10-01: 2,951);
  - IsoCon: 109 outputs, 1,600 reads used (54%);
  - 21 flagged not in the reference; 16 linked back as alleles; 5 new-copy contigs;
  - **3 candidate copies**, 2 counted (>= 2 contigs).
  - The 25 copies all lift (T2d 16, T2i 6, T? 3).
- **Target reads R_T:** 94 net reads fit his maternal chr18 6.72-6.82 Mb best by >= 5 ms. All 29 of register 1360's are among them.

| candidate | contigs | best fit in his genome | class | reads (in R_T) | best reference locus |
|---|---|---|---|---|---|
| iso_NPIP_16 + 17 | 2 | mat chr18 6,766,688-6,803,775, 0.9991 | haplotype-only, **in T** | 8 (7) | NPIP.1 |
| iso_NPIP_14 + 15 | 2 | mat chr18 6,744,425-6,756,028, 0.9986 | unmatched (< 0.999), **in T** | 19 (19) | NPIP.10 |
| iso_NPIP_19 | 1 (not counted) | mat chr18 100.86 Mb, 0.8718 | unmatched | 5 (1) | elsewhere |

- **E1 holds, E2 holds:** exactly one counted haplotype-only candidate lies in T.
- **E3 fails:** it holds 7 of the 94 target reads (bar 47).
- **E4 fails:** 7 of its 8 reads are target reads (0.875; bar 0.90).
- **E5 fails:** the other counted candidate is "unmatched", so it counts as false by the frozen class. But it lies in T as well, 0.0004
  below the 0.999 consensus bar.
- **Where the 94 target reads went:** IsoCon left **66 unassigned (27 of register 1360's 29)**; 26 are in the two target candidates
  (pure: 26 of their 27 reads); 1 is in an output linked back as an NPIP.5 allele; 1 is in the uncounted candidate.

**Reading.** The chain raises no flag outside the maternal-only locus, apart from one single-contig candidate under the floor. It
does find that locus, but splits it in two and recovers only 28% of its reads. The loss is upstream of the flag-link-merge steps:
IsoCon does not converge on most of the locus's reads, as on 10-01 ("real reads too truncated and sparse").

Audit: below (FAIL holds with caveats; the reading above is corrected).

## Audit (metric-auditor, 2026-10-21): FAIL holds with caveats; the reading above is corrected

Every number reproduces:
- **Chain:** net 2,955; 109 outputs; 1,600 reads assigned; 21 flagged / 16 linked / 5 new-copy contigs; 3 candidates.
- **Bars:** E1-E5 = holds, holds, fails, fails, fails.
- **R_T:** 94 reads. The same set comes out at ms margins of 0, 1 and 10; without the intron filter it is 93. Every target read fits T
  best by >= 103 ms.
- **Fates:** 66 / 19 / 7 / 1 / 1. Register 1360's reads: 27 of 29 unassigned.
- **Code path:** only the two declared deviations. The code is byte-identical to the results commit, and the chain is unchanged since
  10-01 apart from docstrings.

**Corrections:**
- **T holds three transcription units, not one copy split in two.**
  - C (6.767-6.804 Mb) is register 1360's copy. It holds 24 of the 29 reads; the other 5 are in A (6.72-6.74 Mb), none in B.
  - B (6.744-6.756 Mb) is a separate unit on the same strand, with its own start site: 21 reads begin at 6,756,009-6,756,064. Its 5′
    exons match RRN3's first three exons (GGO chr18 29.59 Mb). No read joins B and C, and their contigs share no sequence.
  - **Each unit that got a candidate got exactly one:** C → iso_NPIP_16 + 17 (haplotype-only), B → iso_NPIP_14 + 15. A got none.
- **E5's "false" flag is consensus noise on a real unit.** iso_14 and iso_15 each differ from his assembly by ±1 A in poly(A) runs of
  16-17 nt, and the reads covering those runs disagree too (26 / 36 and 22 / 36 carry indels).
- **E1 and E2 hold by one base.**
  - iso_16 scores 1055 / 1056 = 0.99905; iso_17 alone would be unmatched at 0.99898.
  - Had iso_15 matched one base better, E2 would fail instead of E5.
  - The bars could not credit a second true unit in T: a design flaw, since the target window was assumed to hold one copy.
- **E3's bar was out of reach even for a perfect candidate.** No single unit has 47 reads (C 37; B 21, plus 21 that run on into A).
  Per unit:
  - C: 7 / 37 held, 30 unassigned;
  - B: 19 / 21 held;
  - B running into A: 0 / 21;
  - A: 0 / 11.
- **E4 is genuine, but on n = 8.** SRR27438212.5855649 has no record in T and fits pat chr18 22,085,194 at NM 0.
- **Table and text:**
  - "NPIP.10" for iso_14 + 15 came from a 172-kb GGO hit with a 168,474-bp intron; none of its exons touches NPIP.10.
  - "No flag outside T" should read "no new-copy candidate outside T": 16 flagged outputs were linked back.
- **History:** IsoCon converges where a unit's 5′ ends agree (B: 19 / 21) and fails where they spread (C spans 44 kb: 7 / 37). That
  matches 10-01 for C only.

**What stands:** FAIL as registered.
- The chain raised no new-copy candidate outside the maternal-only locus.
- It named the locus's two transcribed units with one candidate each, one of them misclassed by a poly(A) homopolymer error.
- It recovered only 7 of the 37 reads of register 1360's copy.

Detection works here; read recovery depends on IsoCon converging, which fails on transcripts whose 5′ ends spread.
