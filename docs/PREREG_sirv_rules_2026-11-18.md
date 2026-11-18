# Row 1801: four rules from the SIRV diagnosis, tested off the spike-ins (DRAFT of 2026-11-18, bars not yet accepted)

**Status: draft.** The user asked for the two switches (2026-11-18); the bars below are proposed, not accepted, and no arm on a
whole dataset has run. Workspace `r1801_sirv_wrong`.

**Where it comes from.** On the SIRV contigs of the LRGASP mouse ES library (rows 1779 / 1780) ours emits 118 chains, 47 correct
and 71 wrong. Diagnosis of 2026-11-18 on those reads (not registered; after row 1793 the dataset is diagnostic data):
49 wrong chains carry a splice site shifted by 2-6 bp, 16 join true introns in a combination no read carries (15 made by
transfrag extension), 4 are an intron in an unspliced spike-in (2-3 full-length reads with an internal deletion against
282-1,911 reads straight through), 2 other. The default collapse of coordinate duplicates counts 33,954 reads of the correct
chains as 2,885 and 545 reads of wrong chains as 411, so the 0.02 isoform fraction stops removing errors at deep loci. The
shadow rule then drops the true single-exon transcript on the four spike-ins that carry a spliced artefact (4 of 15 lost).

**Rules.**
| arm | change | status |
|---|---|---|
| K | `--keep-coordinate-duplicates` (every read counted) | existing switch; negative as a default on primates (`docs/archive/2026-09/PREREG_primary_dedupe_2026-09-22.md`: gorilla -1,096 matching chains at +3.2 points precision). Re-opened on the user's word |
| E | `RUSTLE_EXTEND_TRANSFRAGS=0` | existing; extension is a default since 2026-11-16 |
| S | `RUSTLE_POLISH_SPLICED_SHARE=0.05` (new, `--polish-spliced-share`) | drop a transcript with an intron spliced by < 5% of the primary reads across it |
| M | `RUSTLE_POLISH_MONO_SHADOW_RATIO=1` (new, `--polish-mono-shadow-ratio`) | the shadow rule keeps a single-exon transcript with at least as many reads as the spliced transcript shadowing it |

On the SIRV contigs (binary `target_r1801`, `arms.sh` / `arms2.sh`; chains / correct / wrong): base 118 / 47 / 71; S 109 / 48 / 61;
K 67 / 50 / 17; K + E 58 / 50 / 8; K + E + S + M 53 / 50 / 3. M alone: chains unchanged, single-exon transcripts on the 15
unspliced SIRVs 11 -> 15, and 4 more on SIRV1-7 (2 in the truth, 2 not). With both switches unset the new binary's `x.gtf`,
`x.final.gtf` and `x.params.tsv` are byte-identical to the old binary's.

**Arms.** base, S, M, K, K + E, K + E + S + M; `assemble` as in row 1793 (`RUSTLE_PTC=0`, same driver mechanism).

**Step 1, shown without a bar** (diagnostic data of row 1793): fly simulated, mouse simulated, mouse ES, OR6737, A119b; SIRV shown.

**Step 2, bars on data not used for a choice:** U1 zebrafish simulated and U2 mouse heart (ENCSR938RPH) of row 1800, which row
1800 scores for other rules but does not mine. Unit: exact intron chains (`bench/famcmp/curve.py`) on `PREFIX.gtf`; `PREFIX.final.gtf`
shown. Per dataset:
- **D1, spliced share.** S against base: precision not lower, and at most 0.5% of base's correct chains lost.
- **D2, shadow support.** M against base: the chains are the same set; of the single-exon transcripts M adds, the share that
  matches a reference single-exon transcript (gffcompare `=`) is at least that share among base's single-exon transcripts.
- **D3, every read counted.** K against base: precision higher, and at most 5% of base's correct chains lost ("pays", row 1793).
- **D4, extension off.** K + E against K: precision higher, and at most 5% of K's correct chains lost.
Gained and lost chains are reported separately for every comparison.

**Decision.** A rule whose bar passes on U1 and U2 is a candidate; the user rules on every default. Nothing changes a default here.

**Limits, stated in advance.** U1 is simulated (no truncated reads, no retained introns: D1 can hardly fail there, D3 and D4 are
flattered). U2 shares the species of one diagnostic dataset and GENCODE is incomplete. The chain score filter was fitted on
de-duplicated counts, so `final.gtf` under K is outside its fit. S reads the BAM once more per contig (wall time reported).
Two primate samples stand for twelve.
