# HANDOFF — state as of 2026-08-15 17:05 PDT

Written after a terminal crash, so that any agent can resume without re-deriving context.
Everything below is either on disk or in git. **Nothing is lost.**

Branch `dna-from-genome`, HEAD **4285b8f**, pushed to origin, **0 unpushed commits**.
Only `tools/stringtie` (third-party submodule) is dirty — leave it.

---

## 1. WHERE THE THREE OBJECTIVES STAND

### O1 — stands. Ship it.
* P1 seed-invariance is a **theorem**; the seed was never in the definition.
* **Both** error rates measured: false-merge **1.33%** [0.37, 4.73] (power measured);
  false-omission **5.6%** [3.0, 10.2] (9/162, unbiased arm; 24.1% pro-edge upper bound).
* Reach **0.5500** [0.3983, 0.6929], chr1 shown representative (Fisher p = 0.6090).
* One live parameter: coverage. τ inert, γ inert on seeded nodes. **Identity has failed 0 times
  across four substrates** (0/728, 194/194, 171/171).
* λ (edge-connectivity) certificate SHIPS in `families.tsv`. Reported, never used to decide
  membership — a 2-copy family has λ=1 necessarily.
* Spec: `docs/seeded_family_definition.md` §1★.

### O2 — the method stands, the OBJECTIVE STATEMENT does not. Fix the prose.
* ⚠ As written ("assign under MAPQ-0 ambiguity") it targets **0.04%** of reads in its own loci and
  duplicates minimap2 on **99.9%** of them. Both measured. Neither survives an examiner.
* What IS defensible, all measured:
  1. the objective **decomposes** ⟹ per-read argmax is the optimum (EM changed **0/3,081**);
  2. **minimap2 is at chance about whether a read belongs at all** — AUC **0.4944**, median MAPQ
     60 vs 60 — while the divergence channel gets **AUC 0.7995**. First non-circular validation;
  3. **12.16%** of tight ties (≤1% AS margin) have a secondary that fits better by divergence.
* ⟹ **Restate O2 as ABSTENTION, not reassignment.** Never claim it assigns better than minimap2.
* ⚠ Still fully open: no ground truth for **reassignment** (the non-circular test covers abstention).

### O3 — a bounded negative, in two strata. Publishable as such.
* Unique/non-collapsible sequence: **M ≤ 6.4** missing expressed copies.
* Collapsible/paralogous: π = 1/35 = 0.0286, 0/26 at cov ≥ 0.8 ⟹ **formally unbounded, vacuous**.
  **O3's target class sits in the dead stratum** — a collapsed copy is ABSORBED, not orphaned.
* One real detection: **STON1 + GTF2A1L**, ~116.7 kb absent from mGorGor1, 125 near-full-length
  unmapped reads, gapless chromosome, 0 GFF lines (COL1A1 control 104), present in chimp+orangutan.
  ⚠ **SINGLE-COPY, n=1, p=0.055 — STRONG BUT UNCONFIRMED. Not a family copy.**
* Docs: `docs/o3_missing_copy_evidence.md`, `docs/o3_unmapped_route.md`.

---

## 2. ⭐ THE BIG UNBLOCK FOUND JUST BEFORE THE CRASH

**The downloaded gorilla DNA is usable whole-genome, NOT Y flow-sorted.**
Artifacts: `/home/juanfra/winloci_scratch/o3_reframe/` (`RESULTS.md`, `PREDICTIONS_prereg.md`,
`kmertest.c`). Predictions were pre-registered before computing.

| sample | R_Y/chr1 | R_X/chr1 | R_chr20/chr1 |
|---|---|---|---|
| SRR26039725/26/27, SRR26152581 | **0.47–0.59** | ~0.50 | ~1.0 |

H0 (male WGS) predicted R_Y ≈ 0.50, band 0.35–0.75 → **all five inside**.
H1 (Y flow-sorted) predicted R_Y ≥ 10 → **rejected by 17× to 1600×**.
Controls valid: R_chr20 ≈ 1.0 (normalisation), R_X ≈ 0.50 (the instrument resolves a 2-fold ploidy
difference, so it would trivially resolve the 26–1600× H1 predicts).

⚠ **Why the metadata lied — do not trust ENA `sample_title` for these runs.** SAMN04003007 is a
2015-era BioSample registered for a Y flow-sorting experiment on KB3781 and later REUSED as the
umbrella BioSample for all mGorGor1 sequencing. SRA auto-composes `experiment_title` from
`sample_title`, so the stale string propagates. The same field labels the **fibroblast IsoSeq** run
`SRR27178662` as "Gorilla Y flow-sorted DNA", which is self-evidently false.

**By-product: `SRR26152581` is HiFi, not CLR** — 8.39% of read 31-mers hit unique genome 31-mers vs
8.42% for Illumina Q30; CLR at ~87% accuracy would give 0.87³¹ ≈ 1.3%, ~6× lower. The `_subreads`
in the filename is a naming artefact; quality strings are Q93-capped.

⟹ **This unblocks two long-open items at once**: the Illumina/HiFi depth check (S2, the field's own
standard, which we never did) and the DNA arm of the reframed O3 below.

---

## 3. WHAT WAS RUNNING WHEN IT CRASHED (both dead, both resumable)

### (a) O1 replication on the matched fibroblast arm — DIED mid-run, must restart
**Why**: the shipped O1 catalog was built on **testis / OR6737** — a different animal from the
reference. The fibroblast arm is **KB3781, the assembly's own cell line**, and O1 has never been
rebuilt on it. Agreement across a different tissue AND individual is a strong robustness claim.

Inputs already prepared (**do not redo**):
* `/mnt/linuxdisk/home/juanfraitu/o1_replicate/fibro_ds.bam` (6.3 GB, indexed) — fibroblast
  depth-matched to testis at **3,962,281 primaries** vs a 3,961,315 target (0.024% off).
* `/mnt/linuxdisk/home/juanfraitu/o1_replicate/build.sh` — the exact invocation.

Progress reached before the crash (from `fibro_gwcat.time`, stderr):
`105,366 skeletons -> 17,564 reps over 26 contigs`; mis-chain filter dropped 247 → **84,015
transcripts** (testis: 289 dropped → 69,528 → 12,415 reps). So at matched depth the fibroblast arm
gives **+20.8% transcripts and +41% reps**. No output files were written.

**Resume**: `cd /mnt/linuxdisk/home/juanfraitu/o1_replicate && nice -n 5 ./build.sh`
⚠ Two bugs in `build.sh` to fix first: the pipeline's progress goes to **stderr**, which the script
sends to `fibro_gwcat.time` (mixed with `time -v` output) while `fibro_gwcat.log` stays empty —
swap the redirects. Runtime was >80 min and unfinished; budget several hours, 4.4 GB steady RSS.

**Then compare on the SHARED NODE SET ONLY** (edge-Jaccard, block agreement, ARI). The two arms will
NOT have the same node set (+41% reps), and comparing partitions over different node sets is
register rule **T7** — the trap that produced a 43% artefact on 2026-08-14.
⚠ Disagreement is **uninterpretable**: tissue and individual are perfectly confounded. Report
"differs, cause unattributable" rather than inventing a story. Agreement is the informative outcome.

### (b) O3 reframing assessment — DIED after phase 1; phase 1 results are on disk
`Workflow({scriptPath: ".../o3-reframe-inter-individual-cn-wf_7948edcd-def.js", resumeFromRunId: "wf_7948edcd-def"})`
Completed agents replay from cache. The Y flow-sorting question (§2 above) was its phase-1 output
and is **already answered** — that was the decisive part.

---

## 4. THE OPEN QUESTION THE USER MUST PUT TO THE ADVISOR

The user suspects they may have misread O3. Two readings:

* **(A) as pursued** — the transcriptome carries copies absent from the genome. Sits in a genuinely
  empty niche (no paper has ever discovered a reference-absent copy from transcriptome data), but
  §1 shows the target class is in the vacuous stratum.
* **(B) the reframing** — between two gorillas, copy NUMBER differs at the genome level (collapse);
  the fibroblast individual matches the reference, the testis one should not. O3 = measure the
  deviation from the reference, especially collapsed copies. (User's analogy: Huntington's.)

**Reading (B) is now EXECUTABLE on one arm** thanks to §2: KB3781 has usable WGS, and the matched
arm is a **null by construction** — any CN difference detected there is a false positive.
⚠ But there is still **no DNA for OR6737**, so the between-individual half cannot be done with data
in hand. And (B) enters a mature field (Bailey, Redon, Sudmant read-depth genotyping, Soto famCN
over SGDP n=269), trading an empty niche for a crowded one.
⚠ The Huntington analogy is probably **misleading** — HD is a within-gene repeat expansion with a
pathogenic threshold; gene-family CN is discrete copies of a whole unit. Closer analogues:
**SMN1/SMN2** (the one classic CNV locus well expressed in fibroblasts), AMY1, C4A/C4B, DEFB4.

---

## 5. NEXT ACTIONS, RANKED

1. ✅**DONE 2026-08-15** — O2 restated as abstention (`copy_assignment_definition.md` §0, both
   objective tables, and the three straggler sites). Commits d53cea2 + 2d2a09f.
1b. ✅**RUN 2026-08-15 — UNDERPOWERED, AND THE ROUTE IS REFUTED** (`docs/o2_reassignment_result.md`).
   n=53 from 2 families (n_eff≈3) vs an expected ~680; filter L1 rejects **29/29** anchors in the
   contested families and **0/53** scored reads are near-ties, so the route cannot reach O2's target
   population. Reassignment on real reads remains UNMEASURED. Superseded design note: — designed and feasibility-measured in
   `docs/o2_reassignment_ground_truth.md`, not yet run. Closes the last Tier-1 gap in O2: structural
   anchors (sequence present in one copy, absent in its sister) label a read's copy of origin WITHOUT
   any scoring decision, so trimming the anchor and re-scoring gives non-circular reassignment
   accuracy on REAL reads. **34/162 families carry an anchor ≥100 bp (median 779 bp)**; median 72
   reads on the lower-depth copy. ⚠ Anchor discovery must use the SHIPPED sensitive tier, not asm20.
2. **Run the Illumina/HiFi depth check (S2)** now that §2 unblocks it — the field's own standard,
   on the matched individual, where excess depth means assembly collapse rather than polymorphism.
3. **Finish the O1 fibroblast replication** (§3a) and compare on the shared node set.
4. **Sweep `c` at 0.45 / 0.55 / 0.75.** `c = 0.50` sits ON its no-change boundary (0.4721, 0.5000],
   not inside a plateau, and the c = 0.75 knee was never run. ⚠ Register rule **T8**: an offline
   re-derivation is a hypothesis generator, NEVER a test — locate candidate values offline, then
   re-run the actual pipeline at those values.
5. **Audit the overlapping-copy cases** — 30/194 two-copy families emit a copy overlapping its own
   sister. Node-definition problem, never audited.

---

## 6. HARD CONSTRAINTS — violating these has cost this project real time

* WSL2, 5 cores / 25 GB. **ONE heavy job at a time.** No nohup-background, no waiter loops.
* **NEVER `pkill -f` / `pgrep -f`** — matches its own shell; has killed 5 sessions. Kill by PID
  after checking `readlink /proc/<pid>/cwd`.
  ⚠ Also note `pgrep -x` fails on `gw_family_catalog`: `/proc/PID/comm` truncates to 15 chars, so
  the process registers as `gw_family_catal`. Monitor by PID, not by name.
* `winloci_scratch` free space is **FICTION** (sparse vhdx on C:). Big outputs → `/mnt/linuxdisk`.
* Build ONLY with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`.
* `cargo test --release --all-targets`, output captured to a FILE (a piped `| tail` returns tail's
  exit code). Baseline **788 passed / 0 failed / 11 ignored**.
* Always `-F 2308` before any per-read statistic.
* **Read `docs/NEGATIVE_RESULTS_REGISTER.md` before proposing an approach** — 628 entries, 307 rated
  high redo-risk, and §2's 20 metric traps are rules, not anecdotes.
