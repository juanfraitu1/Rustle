# The no-deletion control of `o3_candidates` (Amendment 14): C1' FAILS (35/53), C2' HOLDS (0.54%) — the default flip 1f49d0f0 is reverted (d04b6ae9) — 2026-10-03

Prereg: `docs/PREREG_rna_allele_haplotype_count_2026-10-01.md`, Amendment 14 (`750f59ee`, written before the run), with Amendment 9's
procedure (classification of each candidate against KB3781's own haplotypes, arm C). Rulings R22 (the default flip 1f49d0f0 ships iff
C1' and C2' hold, else it is reverted) and R23 (the stage's cost on a full BAM is measured here). Stage: `o3_candidates` rebuilt from
HEAD `75827d5d` (sha1 `77eed47a...`; its source differs from the A13 binary's `0f5824a7` in comments only). Recipe:
`bench/rna_allele/accept_o3_candidates.sh` with `ACC=a14` (+ `.py`, `panel_to_copies.py --all`; commits `b272bc46` + `50cf5864`, the tie breakdown and the a14 guard). Work dir
`/mnt/linuxdisk/tmp/rna_allele/a14/` (the whole-BAM run in `a14/wholebam/`). Scorer and helper outputs copied to
`docs/O3_CANDIDATES_CONTROL_A14_score.out.txt`. The run under control: `docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md` (A13); the
IsoCon chain's control: `docs/RNA_ALLELE_CONTROL_2026-10-01.md` (A9).

## Verdicts

| rule | registered bar | measured | verdict |
|---|---|---|---|
| **C1'** (specificity at the stage's operating point) | families with >= 1 false flag (b allele, c unmatched; `pri` counts too) <= 1/3 of A13's detection 25/53, i.e. **<= 8 of 53** | **35/53 = 66.0%**; family-level LR of a flag (25/53) / (35/53) = **0.71** | **FAILS** |
| **C2'** (cost without a deletion) | reads placed on a candidate not derived from their own copy <= 5% of all reads | **318 / 59,013 = 0.54%** | **HOLDS** |
| decision (Amendment 14, R22) | the flip ships iff C1' and C2' hold | C1' fails | **1f49d0f0 REVERTED in d04b6ae9: `candidates` is opt-in again (R14)** |

Reported beside (Amendment 14): families with any flag (a + b + c) **35/53 = 66.0%** (the same 35: the one family with a true flag,
GWFAM175, also carries false ones); A9's IsoCon chain at any support: **16/53 = 30.2%** with a false candidate.

## Substrate and procedure (Amendment 9's, the stage in place of the IsoCon chain)

- The 53 families of Amendment 7, all **201 copies, nothing masked**; the 59,013 scored reads as `control/R0.bam` (Amendment 9's arm R0:
  aligned to the unmasked `_pri` with the pipeline's flags; no unmapped record); `GGO.fasta` + `winloci_data/GGO.splice.mmi`.
- Copies table `A14.copies.{tsv,fa}`: `panel_to_copies.py --all` over `mask` + `keep` of `linktest/panel.json` (the masked copy first,
  Amendment 9's order), sequences from the unmasked `_pri`, `n_reads` by `R0.bam` (every copy has a primary read).
- The stage at its defaults (`--min-support 6`, delta 0.00958, `--min-cluster 3`, the 1,000-read cap, `--threads 4`, no
  `RUSTLE_CACHE_DIR`), batched as A13 (A12's five batches, `a12/batches.txt`).
- **Classification (C1').** The 56 flagged unions, renamed `iso_<family>_<k>`, aligned to KB3781's maternal and paternal splice indexes
  with `minimap2 -c -x splice:hq -uf -N 20` (Amendment 9's command); each takes the class of its best hit by identity x coverage
  (`control_test.py classify`'s rule, mirrored in `accept_o3_candidates.py classify`): **a** = >= 0.999 on the haplotype `_pri` did not
  take that chromosome from (`chrmap.tsv`), outside the lifted B interval of every copy of the family (`control/copies_lift.tsv`,
  lift_frac >= 0.5) — a true flag; **b** = inside one — an allele; **c** = no hit >= 0.999 — unmatched; **pri** = >= 0.999 on `_pri`'s own
  haplotype. Self-check: the harness re-classifies Amendment 9's 28 candidates from Amendment 9's own PAFs and reproduces
  `control/candidates.tsv` exactly (24 c, 3 b, 1 a).
- **Labels ("derived from copy g").** Each union's best hit in the unmasked `_pri` (`-c -x splice:hq -uf -N 20`) overlapping copy g's
  clean interval (A9's copy order): 55 of the 56 unions carry a copy of their own family, 1 `elsewhere`.
- **Arm C (C2').** `_pri` + the 56 unions (`C.fa`, `minimap2 -x splice -d`); the three scored parts realigned with R0's flags
  (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, Amendment 9's arm C command) and merged; each union its own locus (the
  stage's candidate is already its component); reads classified as `control_test.py score` does (a 0.98 AS tie across loci is
  unplaced). Self-check: the harness's arm R0 equals Amendment 9's `control/score.json` exactly.

## The stage without a deletion

| batch | families | `/usr/bin/time` Elapsed | stage clock (nets / clusters / genome) | peak RSS | flagged (families) |
|---|---|---|---|---|---|
| 0 | 6 | 171.6 s | 153.7 s (4.7 / 89.2 / 59.8) | 15.7 GB | 3 (2) |
| 1 | 12 | 218.8 s | 197.3 s (3.9 / 149.4 / 43.8) | 15.7 GB | 12 (7) |
| 2 | 12 | 394.9 s | 355.3 s (4.3 / 242.6 / 106.6) | 15.3 GB | 12 (9) |
| 3 | 11 | 302.2 s | 273.5 s (4.3 / 193.6 / 75.1) | 15.9 GB | 15 (9) |
| 4 | 12 | 256.5 s | 231.3 s (3.8 / 176.3 / 49.8) | 15.1 GB | 14 (8) |
| **sum** | 53 | **1,344.1 s = 22.4 min** | 1,211.1 s | 15.9 GB (max) | **56 (35)** |

- **Products:** nets 59,150 reads (44,782 after the cap); 414 clusters — 175 in the reference, 159 linked, 80 new copy — -> 79
  candidates, **56 flagged (>= 6 reads) in 35 families** (A13, one copy per family masked: 109 / 82 in 41). Clusters per flagged
  candidate: 1 in 55, 2 in 1 (`iso_GWFAM244_0`).
- **Phase-1 counters (stage logs, 53 families):** 307 read clusters >= `--min-cluster`; 0 empty consensus dropped; refinement split off
  3,162 reads into 132 new clusters (0 fell under `--min-cluster`); 13 kept sets re-templated (12 families: GWFAM28, 149, 182, 247, 269,
  272, 317, 331, 348, 398, 402, 490); the significance merge absorbed 25 clusters in 70 rounds, 0 absorptions undone; fallback
  templates (no eligible member): 16 the longest member with an aligned partner, 4 the longest member; 414 final clusters.
- **Attribution (pass B) without a deletion** (stage logs; `nets` re-runs it per batch and reproduces every count and n_net / n_used of
  all 53 families): unmapped >= 300 bp **0** in every batch (R0 has no unmapped record); poorly placed >= 300 bp 1,314 / 461 / 1,269 /
  1,059 / 1,289 = 5,392 (8 below the floor); aligned 4,776; attributed 1,348, of which 1,303 to a family of ANOTHER batch (1,301 =
  99.8% their own family: reads their family's batch nets by pass A, ruling R18). **Joined a net of their batch: 45 joins of 37 reads,
  all to a wrong family** (GWFAM169 20, GWFAM401 17, GWFAM244 7, GWFAM164 1) — by construction: a read of a family in the batch is in
  that family's net already, so with nothing missing only another batch's reads are eligible. **None reaches a flagged candidate's
  cluster** (16 in linked clusters, 24 in no reported cluster, 5 over the cap; in arm C 42 are unplaced and 3 stay on their own copy):
  the false flags do not come from the attribution.

## C1': the flags against the diploid truth — FAILS

| class | flags | families |
|---|---|---|
| a, haplotype-only locus (true) | 1 | GWFAM175 |
| b, allele | 1 | GWFAM175 |
| c, unmatched | 54 | 35 (all) |
| pri | 0 | – |

- **The one true flag is Amendment 9's reference-absent copy:** `iso_GWFAM175_1` (183 reads) at identity x coverage 1.0000 to paternal
  chr5 `CM054563.2:40,027,838-40,031,943` (A9's 12-transcript candidate: `40,028,172-40,031,100`; register row 1204), outside the lift of
  every copy of the array. Its siblings are false: `iso_GWFAM175_2` (116 reads) is an allele (1.0000 at paternal
  `CM054563.2:41,468,509-41,471,848`, inside the lift of GWFAM175:5) and `iso_GWFAM175_0` (191 reads) is unmatched.
- **The 54 unmatched flags exist in neither of KB3781's haplotypes at 0.999.** Over these 54: best haplotype hit median 0.956
  (0.486-0.9985; 2 >= 0.99, 30 >= 0.95, 11 < 0.90); whole-length d to the nearest `_pri` locus median 0.045 (0.0097-0.515; 13 <= 0.02,
  12 > 0.10); reads per flag median 196.5 (6-534; 39 with >= 100 reads, 8 with <= 10) — well supported, not few-read noise. (Over all
  55 false flags, b + c: reads median 191, 40 with >= 100; d median 0.044.) What they are — mosaic consensus of paralogs, transcript
  structures absent from both assemblies, or residual consensus error — is not determined here.
- **They recur:** 28 of the 56 flags reproduce one of the A13 run's 46 survivor-derived flags (same family, identity >= 0.999 over >= 90%
  of the shorter sequence: Amendment 9's overlap rule) — 28 of those 46 arise with nothing deleted — and 29 of the 31 families with an
  A13 survivor-derived flag carry a false flag here; 18 of the 25 families whose deleted copy A13 found do too. Against A9's 28 IsoCon
  candidates: 7 flags reproduce one (7/28 reproduced); 10 families are flagged by both (A9 16, A14 35).
- **C1': 35/53 = 66.0% > 15.7% (8 families) -> FAILS.** The family-level likelihood ratio of a flag is (25/53) / (35/53) = 0.71: at the
  stage's operating point (`--min-support 6`) a flag in a family is no evidence that one of its copies is missing. (A9's IsoCon chain
  failed its C1 by two families, 16/53 against 14.)

## C2': arm C — HOLDS

| | R0 (`_pri`) | C (`_pri` + the 56 flagged unions) | A9's C (`_pri` + IsoCon's 28 candidates) |
|---|---|---|---|
| stay on own copy | 58,097 | 50,146 | 55,306 |
| on a candidate derived from own copy | – | 3,774 (6.4%) | 2,018 (3.4%) |
| **false move** (a candidate not derived from own copy) | – | **318 (0.54%)** | 27 (0.05%) |
| other copy of the family | 13 | 12 | 11 |
| elsewhere | 0 | 0 | 0 |
| unplaced | 903 | 4,763 | 1,651 |

- **C2': 318 / 59,013 = 0.54% <= 5% -> HOLDS.** 295 of the 318 are GWFAM175 reads on `iso_GWFAM175_1` — the real paternal-only copy
  (class a), whose reads the labels can only give to a `_pri` copy (as in A12's and A13's arm M) — then 19 GWFAM244 reads on its
  `elsewhere`-labelled `iso_GWFAM244_1`, 2 GWFAM175 reads on `iso_GWFAM175_2`, 1 + 1 on GWFAM331's and GWFAM37's flags.
- Beside (no rule): 4,165 reads placed in R0 become unplaced in C (by family: GWFAM47 707, GWFAM98 476, GWFAM173 446, GWFAM158 433,
  GWFAM54 410, GWFAM269 375, GWFAM348 342, GWFAM268 214, ...). By their 0.98 AS ties in C: 4,074 tie between their own copy and a false
  flag derived from it; 81 tie their own copy with a flag derived from another copy (or labelled `elsewhere`; 43 of them with GWFAM402's flag);
  10 do not tie with their own copy (9 tie between two flags). 3,774 move onto a union derived from their own copy; stay falls 13.7%.
  In O2 every one of the 55 false flags is an extra copy.

## The decision (Amendment 14, ruling R22): the flip is reverted

C1' fails, so the default-on flip (1f49d0f0) does not ship. **Commit `d04b6ae9` reverts it** with the A13 doc's revert recipe (three
files, five hunks: `tools/rustle_pipeline.sh` takes the pre-flip comment above the `case "$STAGE"` block; `README.md` the pre-flip
pipeline paragraph; the spec's header, §9 and §9b keep the record of the flip and R22-R24 with a dated outcome line), the outcome and the
measured R23 cost written where the recipe re-added "not yet measured" (the README's cost clause and the driver header cite "spec §9b,
R23"). Everything else 1f49d0f0 changed reverted without conflict and equals `aad7aca3`. Checks: the pre-flip driver e2e
(`tests/fixtures/o3_candidates/driver/run_e2e.sh`, opt-in check (d)) passes 9/9; `cargo test --release --lib module_status` 2/2. The
stage stays available opt-in (`--candidates`, ruling R14). A13's verdicts (A13-1/2/3) stand unchanged; the control's numbers are the
next prereg's starting point.

## The stage's cost on a full BAM (ruling R23): one batch of 50 families did not finish in a 10-minute call

**Setup.** The full gorilla fibroblast Iso-Seq BAM (`fibroblasts/GCA_029281585.2_flnc_mm.bam`: 23.2 GB, 34,860,395 mapped records + 959
unmapped; minimap2 2.31 `splice:hq` to `_pri`); copies table = the 2026-08-14 interval table, 378 families / 915 copies
(`refabsent/panel.json` via `panel_to_copies.py --all`, sequences from `_pri`; every copy has primary reads, 1,570,327 primary-read
overlaps in all) in `wholebam/W.copies.*`; 8 batches of <= 50 families in table order (`wholebam/batches.txt`); `--threads 4`, the
defaults, no cache; `GGO.splice.mmi` for the genome phase. Batch 0 (GWFAM0 ... GWFAM149 in string order, 143 copies, 630,418 primary-read overlaps; the
other seven batches have 61-226 k) ran twice and was stopped both times; the measurement stopped there, as R23's task directed.

| | attempt 1 | attempt 2 |
|---|---|---|
| stop | the lock's `timeout 585` (SIGTERM; `/usr/bin/time` dies with the stage: no stats) | `timeout -s INT 570` (`/usr/bin/time` survives SIGINT and reports) |
| elapsed at the stop (realtime) | 650 s (bash `time` 649.8 s; terminal output, saved in the score copy) | **631.2 s** (`/usr/bin/time` Elapsed) |
| pass A + the pass-B sweep done (`attrib.fa` complete) | +155 s | +146 s |
| attribution targets written (`attrib_targets.fa`) | +183 s | +218 s |
| attribution alignment done (`attrib.paf` mtime) | +396 s | +468 s |
| the stage's pass-A / pass-B log lines = the end of the nets phase (log timestamps) | **+397 s** | **+469 s** |
| families clustered when stopped (of 50) | 14 (16.4 s each) | 7 (19.9 s each) |
| peak RSS (the stage and the minimap2 runs it had reaped) | – | **10.7 GB** (10,688,296 kB); user 1,557.7 s, sys 164.8 s |

- **Pass A:** 627,388 reads by a primary record and 418,841 secondary records on the batch's copies; pass B found the primaries of
  20,356 secondary-only reads.
- **The attribution set** (written during the sweep): **88,571 reads, 167 MB** — 202 unmapped >= 300 bp + 88,369 poorly placed >= 300 bp
  in no net of the batch (206 below the floor). **The targets: 649,861 records, 2.61 GB** — the batch's 648,946 net reads (a read in two
  nets once per net) + the 915 copies (record counts from the terminal output of attempt 1's temp dir before it was deleted, saved in
  the score copy; the size from its listing). Aligned 61,016 (117 unmapped, 60,899 poorly placed); attributed 30,457 (5 / 30,452); **joined the
  batch's families 26,710** (no truth labels on this library: right / wrong unknown).
- No flagged candidate: the genome phase was never reached.
- **Cost statement.** On a full library one call of 50 families spends 397-469 s on its nets alone (start to the pass-A/B log lines)
  — a sequential sweep of the whole 23-GB BAM plus one alignment of ~89 k reads against the batch's whole nets — before any
  clustering; at 16.4-19.9 s per family the batch would need about 20-25 min (397-469 s + 50 x 16.4-19.9 s = 20-24 min before the
  genome phase: an estimate from the measured phases, not a measurement), and every batch repeats the sweep. The stage
  does not fit 10-minute calls at ~50 families per batch, and its nets phase scales with the library (the attribution set) and with the
  batch's expression (the targets).
- Clock note: the realtime clock (`date`, gawk's `systime()`, `/usr/bin/time`) ran ~11% ahead of the timeout's timer (585 s -> 650 s, 570
  s -> 631 s), as A12's and A13's docs noted for the stage's own clock; the verdict ("did not finish in 10 minutes") holds on either.

## Post hoc (not the verdict; starting points for the next prereg)

- **A read floor does not separate:** families with a false flag at >= 10 / 20 / 50 / 100 reads per flag: 32 / 31 / 30 / 28.
- **The >= 2-cluster floor** (ruling R1's alternative, reported beside in A12 and A13) leaves 1 false-flag family (GWFAM244) against
  A13's 7/53 detection at that floor (LR 7) — but that floor finds only 7 of the 53 deleted copies (A13: 25 at `--min-support 6`).
- The survivor-derived flags A13 reported as a cost beside its verdict (12 -> 46) are mostly these same locus-specific false flags (28 of
  46 recur with nothing deleted): the next rule has to separate a union of a copy's own reads that misses both haplotypes from a copy
  that is not in the reference.

## Caveats

- One individual (KB3781, the reference animal), fibroblast Iso-Seq, the 1,000-read cap. The a / b split rests on Amendment 9's asm5 lift
  (200 of 201 copies lifted; GWFAM175:2 did not) and, as in A9, on single PAF records at 0.999 (a union split across records counts as
  unmatched).
- Batching (ruling R18) changes only the attribution here (45 wrong-family joins, none in a flagged candidate's cluster); clusters,
  candidates and genome hits are per family.
- The binary under test was rebuilt from 75827d5d (sha1 `77eed47a...`) because the one at the path (sha1 `6775fddb...`, rebuilt at 21:18
  by the Task 3 fix round's `cargo test --release`) was not the A13 binary (`0c97f623...`); `git log 0f5824a7..HEAD -- src/` shows
  comment changes only, so the stage's logic is A13's.
- The whole-BAM figures come from one batch, stopped twice; the per-phase times varied between the attempts (targets 28 vs 72 s,
  alignment 213 vs 250 s); the full eight-batch cost was not run.

## Appendix: per family (the 35 flagged families)

Flags as class (a / b / c), reads, whole-length d, best haplotype hit (identity x coverage). "A9 / A13-S": flags reproducing one of A9's
candidates / one of A13's survivor-derived flags. Arm C counts are of the family's reads.

| family | batch | net (used) | clusters: ref / linked / new | flagged: class, reads, d, best haplotype hit | A9 / A13-S match | arm C: own-candidate / false moves / newly unplaced | A13: deleted copy found; survivor-derived flags |
|---|---|---|---|---|---|---|---|
| GWFAM4 | 0 | 400 (400) | 0 / 2 / 1 | c 321 0.106 0.8939 | 0 / 1 | 0 / 0 / 10 | no; 1 |
| GWFAM6 | 1 | 1060 (1000) | 3 / 2 / 3 | c 457 0.170 0.8296 | 0 / 0 | 81 / 0 / 6 | no; 2 |
| GWFAM23 | 3 | 2000 (1000) | 3 / 5 / 1 | c 234 0.083 0.9183 | 0 / 0 | 241 / 0 / 0 | yes; 2 |
| GWFAM37 | 4 | 1315 (1000) | 0 / 6 / 2 | c 369 0.035 0.9649; c 361 0.515 0.4856 | 0 / 0 | 85 / 1 / 20 | yes; 1 |
| GWFAM47 | 1 | 2048 (1000) | 7 / 3 / 4 | c 224 0.043 0.9571; c 164 0.090 0.9108; c 94 0.024 0.9759 | 0 / 0 | 103 / 0 / 707 | yes; 4 |
| GWFAM54 | 2 | 1082 (1000) | 4 / 0 / 1 | c 443 0.081 0.9191 | 0 / 1 | 30 / 0 / 410 | no; 1 |
| GWFAM98 | 4 | 1500 (1000) | 5 / 1 / 2 | c 7 0.015 0.9823 | 1 / 1 | 13 / 0 / 476 | yes; 1 |
| GWFAM99 | 1 | 1500 (1000) | 2 / 2 / 3 | c 317 0.055 0.9447; c 309 0.408 0.6375 | 0 / 0 | 212 / 0 / 9 | no; 2 |
| GWFAM100 | 0 | 1323 (1000) | 3 / 3 / 3 | c 324 0.038 0.9614; c 233 0.040 0.9602 | 0 / 0 | 310 / 0 / 17 | yes; 0 |
| GWFAM104 | 2 | 1312 (1000) | 3 / 2 / 2 | c 534 0.214 0.7847; c 342 0.036 0.9641 | 0 / 1 | 103 / 0 / 68 | yes; 1 |
| GWFAM125 | 2 | 990 (990) | 4 / 1 / 1 | c 449 0.016 0.9845 | 0 / 1 | 166 / 0 / 22 | no; 1 |
| GWFAM158 | 4 | 1726 (1000) | 8 / 8 / 2 | c 223 0.083 0.9182; c 42 0.052 0.9478 | 0 / 0 | 1 / 0 / 433 | yes; 2 |
| GWFAM161 | 4 | 346 (346) | 6 / 5 / 1 | c 14 0.014 0.9950 | 0 / 0 | 8 / 0 / 6 | no; 0 |
| GWFAM163 | 3 | 1187 (1000) | 1 / 8 / 3 | c 9 0.016 0.9761; c 6 0.044 0.9486 | 2 / 2 | 23 / 0 / 4 | no; 2 |
| GWFAM169 | 1 | 572 (572) | 4 / 1 / 1 | c 123 0.016 0.9840 | 0 / 0 | 3 / 0 / 29 | no; 0 |
| GWFAM173 | 3 | 1085 (1000) | 9 / 4 / 3 | c 247 0.100 0.9412; c 191 0.014 0.9856; c 130 0.056 0.9439 | 0 / 1 | 32 / 0 / 446 | no; 2 |
| GWFAM175 | 1 | 2463 (1000) | 5 / 1 / 7 | c 191 0.046 0.9538; **a 183 0.083 1.0000**; b 116 0.023 1.0000 | 2 / 2 | 517 / 297 / 191 | yes; 2 |
| GWFAM181 | 2 | 966 (966) | 2 / 3 / 1 | c 489 0.378 0.6215 | 0 / 1 | 224 / 0 / 13 | no; 1 |
| GWFAM185 | 1 | 415 (415) | 3 / 2 / 1 | c 185 0.010 0.9898 | 0 / 1 | 3 / 0 / 5 | no; 1 |
| GWFAM236 | 2 | 1500 (1000) | 4 / 2 / 1 | c 100 0.044 0.9558 | 0 / 0 | 47 / 0 / 0 | no; 0 |
| GWFAM244 | 3 | 1706 (1000) | 2 / 5 / 4 | c 71 0.068 0.9127 (2 clusters); c 7 0.061 0.9383 | 1 / 0 | 125 / 19 / 1 | yes; 1 |
| GWFAM246 | 4 | 1500 (1000) | 3 / 1 / 3 | c 231 0.023 0.9770 | 0 / 0 | 38 / 0 / 3 | no; 1 |
| GWFAM268 | 4 | 795 (795) | 0 / 3 / 3 | c 228 0.194 0.8060 | 0 / 1 | 3 / 0 / 214 | yes; 1 |
| GWFAM269 | 2 | 1286 (1000) | 3 / 4 / 2 | c 60 0.021 0.9843; c 6 0.010 0.9884 | 0 / 0 | 210 / 0 / 375 | yes; 0 |
| GWFAM282 | 3 | 1500 (1000) | 3 / 1 / 1 | c 312 0.027 0.9729 | 0 / 1 | 227 / 0 / 0 | no; 1 |
| GWFAM317 | 3 | 558 (558) | 7 / 4 / 2 | c 158 0.028 0.9590; c 6 0.014 0.9985 | 0 / 2 | 52 / 0 / 84 | yes; 2 |
| GWFAM331 | 2 | 430 (430) | 1 / 2 / 2 | c 134 0.065 0.9729; c 123 0.122 0.9016 | 0 / 2 | 109 / 1 / 30 | yes; 2 |
| GWFAM335 | 3 | 765 (765) | 3 / 3 / 3 | c 293 0.162 0.8378; c 9 0.187 0.8135 | 1 / 2 | 266 / 0 / 5 | no; 2 |
| GWFAM348 | 4 | 921 (921) | 6 / 3 / 4 | c 338 0.015 0.9888; c 205 0.053 0.9522; c 183 0.012 0.9853; c 59 0.039 0.9654 | 0 / 3 | 276 / 0 / 342 | no; 3 |
| GWFAM402 | 2 | 527 (527) | 4 / 2 / 1 | c 202 0.228 0.7716 | 0 / 1 | 27 / 0 / 51 | yes; 1 |
| GWFAM407 | 3 | 1467 (1000) | 6 / 4 / 2 | c 325 0.012 0.9882 | 0 / 0 | 98 / 0 / 3 | yes; 0 |
| GWFAM425 | 4 | 1082 (1000) | 1 / 3 / 2 | c 349 0.506 0.4944; c 120 0.072 0.9276 | 0 / 1 | 29 / 0 / 18 | yes; 1 |
| GWFAM439 | 3 | 511 (511) | 3 / 8 / 1 | c 8 0.049 0.9555 | 0 / 1 | 2 / 0 / 1 | no; 1 |
| GWFAM440 | 2 | 257 (257) | 3 / 1 / 1 | c 32 0.010 0.9894 | 0 / 1 | 0 / 0 / 21 | yes; 1 |
| GWFAM490 | 1 | 917 (917) | 1 / 5 / 2 | c 262 0.035 0.9671 | 0 / 1 | 110 / 0 / 145 | yes; 1 |

No flag (18): GWFAM21, 28, 62, 105, 112, 123, 144, 149, 164, 177, 182, 192, 247, 272, 314, 398, 401, 415.

Register rows 1226-1230.
