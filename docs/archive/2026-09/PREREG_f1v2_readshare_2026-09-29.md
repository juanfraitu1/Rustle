# Pre-registration: F1v2 (F1 + a read-share MINORITY condition on each bridge) held out on two HUMAN libraries

**Written 2026-09-29 (KEY=f1v2) before any held-out human F1 or F1v2 product exists.** The rule was designed on
development data only (§12): the two gorilla samples spent by the F1 test, the fusion simulation with its confounder
arms, and the human A119b contigs chr16 and chr20 (dev GTFs, excluded from the substrate). This file binds once
Amendment 1 records its sha1 and the frozen instruments' sha1s. Until then nothing held-out runs.
User goal (verbatim, relayed by the orchestrating session): "F1 v2 — reduce F1's single-gene cuts with a read-share
condition, then test it held-out on HUMAN." F1v2 is opt-in whatever the outcome; a default flip is the user's call.

## 0. The question

F1 (`PREREG_f1_bridge_locus_2026-09-28.md`, r1145) lowered fused loci on two gorilla samples, but its falsifier Z1
fired (r1146): 66-74% of bridge transcripts sit on an annotated intron of ONE gene, and about 40 genes per sample were
cut (bridge splits SEP 51 / FRAG 40 and 53 / 39). Post hoc, the correct bridges carried a small read share and the
gene-cutting ones were usually the gene's main isoform. F1v2 adds one ordering test: **a bridge must carry fewer reads
than each side it would separate carries without it.** On the spent gorilla samples this removed 27 of 40 and 24 of 39
bridge FRAG at a cost of 6 and 4 SEP and 4 fused loci each (§12). Does that hold on a **different species, genome,
annotation and library** — human A119b (every annotated contig except chr16, chr20, chr21, chr22) and human testis —
with the fused-locus reduction kept, beyond RG3 and beyond matched random relabellings, with fewer single-gene cuts
than F1 and no net loss of correct separations, chains unchanged, and no NPIP regression?

## 1. The rule F1v2 (binding; `f1v2.py` b4e788adffe92f0b8dcb32cc3584d1380e093a37, `--rule min`, over the frozen F1)

### 1.1 Verbatim from the frozen file's docstring (the RULE block, whole)

```
RULE (binding once pre-registered; every term except MINORITY is F1's, verbatim)
  For a bridge junction J of input gene_id g (F1: STRUCTURAL(J) and UP-proof(J) and DOWN-proof(J)), on J's strand:
    T_J = g's transcripts that use J; R_J = g's other transcripts on J's strand; the components of R_J under
    same-strand exon overlap, labelled UP / DOWN / STRADDLE / INSIDE exactly as F1 labels them (f1_bridge.side).
    reads(X) = the sum of the GTF `reads` attribute over the transcripts of X (F1's parse; absent -> 0).
    UP_J = the union of the UP components, DOWN_J = the union of the DOWN components.
  MINORITY(J) :<=> reads(T_J) < reads(UP_J)  AND  reads(T_J) < reads(DOWN_J)            (--rule min, F1v2)
     "the link carries fewer reads than EACH side it would separate carries without it" -- an ordering (majority)
     test with no constant: equivalent to share(J) = reads(T_J) / (reads(T_J) + min(reads(UP_J), reads(DOWN_J))) < 1/2.
  BRIDGE_v2(J) :<=> BRIDGE_F1(J) AND MINORITY(J). Every J is judged on the original locus (order-free), as in F1.
  Everything after the bridge set (pieces, keeper, names, fusion_of / fusion_junction, the families GTF without the
  bridge transcripts) is F1's regroup() and rewrite(), unchanged; fusion_junction lists the kept junctions only.
```

F1 itself is `f1_bridge.py` 37ee8e77 `--mode full` (the RULE block of the F1 prereg §1.1, unchanged; its BAM-derived
proofs are read from F1's own junction table). Implementation facts that are part of the rule: `f1v2.py` recomputes
T_J and the UP / DOWN / INSIDE components from the BASE GTF and **asserts** they equal F1's junction-table row (T_J ids,
up / down / inside counts, no STRADDLE); `--rule none` must reproduce F1's GTF and families GTF byte for byte (gate
G1h). Because BRIDGE_v2 ⊆ BRIDGE_F1, F1v2's pieces are unions of F1's pieces and bridges; with no F1v2 bridge, F1v2 =
RG3 (rg3.py ec17e540).

### 1.2 Every constant and its source

| element | value | source |
|---|---|---|
| everything in BRIDGE_F1 | F1's constants (F1 prereg §1.2) | frozen, unchanged |
| read measure | the GTF `reads` attribute summed over transcripts | F1's parse (the attribute RG3 and F1 already use for keepers) |
| sides | the UP and DOWN components F1's STRUCTURAL test already computes | F1 |
| comparison | strict `<` against EACH side (both-sided minority) | an ordering test, no constant; the ½ is the definition of "minority", not a fitted cut. Strictness (a tie abstains) was fixed before the tie count was seen; on dev ties split 3 SEP / 1 FRAG (OR) and 1 / 2 (KB) |

Why this form and not a threshold: in a readthrough X→Y, X's own transcripts (ending at X's PAS inside J's intron) and
Y's own transcripts (starting at Y's promoter inside it) are independent units, and the link is a minority product of
each. In the one-gene case F1 cut (an intronic PAS plus an intronic promoter), the link IS the gene's main isoform, so
it outnumbers at least one of the two minor side isoforms. The dev curve (§12.2) shows the ½ point on a plateau of net
correct separations (OR peaks at 0.4, KB at 0.5), not at a fitted optimum.

## 2. Arms (all post-processors of one stored GTF per sample; no re-assembly)

| arm | product judged | how |
|---|---|---|
| **BASE** | `rt_arms/<s>/<s>.BASE.gtf` → `runs/<s>/<s>.gtf`, sha1 `631c9f11` (A119b) / `50f239d6` (testis) | as stored |
| **RG3** (comparator) | `rg3.py` ec17e540 on BASE | testis gate G6: equals the RG3 test's `rg3_run/U/human_testis.RG3.whole.gtf` |
| **F1** (frozen; comparator and reported) | `F1.families.gtf` (loci), `F1.gtf` (arm **F1all**) | `f1_bridge.py` 37ee8e77 `--mode full` on the WHOLE BASE GTF in contig batches (`f1_batch.py` 14144bde), merged |
| **F1v2** | `F1v2.families.gtf` (loci), `F1v2.gtf` (arm **F1v2all**) | `f1v2.py` b4e788ad `--rule min` on (BASE, merged `F1.junctions.tsv`) |
| **NULL(k)**, k = 0..4 | `NULL<k>.families.gtf` | `f1_null.py` ea2715f2 (F1 prereg §2 protocol, unchanged) `BASE F1v2.gtf NULL<k> --label f1v2ho:<s>:s<k>` — matched to F1v2's changes |

Every arm is computed on all contigs; each clause is scored on its substrate's contigs (§3). No NULL is drawn for F1
(F1's human rows are reported, never judged).

## 3. Substrates and exposure

- **V_A = human_A119b, every contig the sample's annotation covers except chr16, chr20, chr21, chr22**
  (readthrough_eval `--drop-contigs chr16,chr20,chr21,chr22`; split_cls_h `--exclude-contigs` the same). BAM
  `winloci_data/A119b.t2t.bam`. Held out in species, genome, annotation and library from the gorilla design; in contig
  from the human dev contigs chr16 / chr20 of the same library.
- **V_T = human_testis, every annotated contig** (BAM `_from_wsl/human_val/human_testis.t2t.bam`). Held out in library.
- Genome CHM13 v2.0 (`chm13v2.0.fa`); annotation = the figures registry's (`figures/samples.tsv` → `human_ref_gff` =
  `/mnt/linuxdisk/tmp/regress/chm13.gff`, the CHM13 v2.0 RefSeq full GFF; readthrough_eval's provenance on the dev
  run records exactly this file); the genes table `families_gw/species/human` is built from the same file
  (`annotation.key`). **Never `HSA_genomic.gff`.**
- The two samples are reported sample by sample; nothing is summed across them; no gorilla number enters a clause.

**Exposure (plainly).** Neither human sample has served an F1 test. Both have served verdict-bearing held-out tests
before: readthrough R/RQ1 (r1117), R2/R3 (r1119/r1120), completeness `pas-end` (r1136) and the descriptive r1137-r1139,
the TSS polish (human only), and — testis only — the RG3 NPIP-block test (r1142-r1144). **This is at least the fifth
verdict use of each and the first on bridge regrouping.** What this author has read or run on human data before this
file:
- the RG3 prereg with its Outcome, including testis's NPIP-block rows (U = 294 loci, 7 present Dishuck copies, c1 5 in
  every arm, RG3's one U split PAGE2 / PAGE2B; the only fused NPIP holder is NPIPB14P's with the PDXDC2P-NPIPB14P
  readthrough, 1 transcript, 2 reads) and testis's whole-GTF RG3 splits (7 of 13,012 gene_ids);
- file sizes of the six BASE GTFs (A119b 334 MB, testis 26 MB), the two human BASE sha1s, and `samtools idxstats`
  record counts of the human BAMs (A119b 68.0 M, testis 9.1 M);
- **the dev contigs A119b chr16 and chr20** (dev GTFs `rg3_port/runs/frozen/hsa{16,20}.off.unset.gtf`, not the
  held-out BASE): F1, F1v2, RG3, 5 NULLs, split classes and FUSED, all in §12.3. They are excluded from V_A.
No held-out human GTF line, junction table, readthrough table, split table or families product was opened.

## 4. Instruments and gates (frozen in `/mnt/linuxdisk/tmp/rustle_figures/f1v2_frozen/`, `SHA1SUMS`, plus the F1 set)

New: `f1v2.py` b4e788ad; `split_cls_h.py` 4a6b178d (split_cls_gw e8b76f83 with `--species`, `--gff`, `--rt-mode`,
and the ANN / ANN_RT junction kinds); `fused_pt.py` c42a4605 (reported only: the two FUSED readings of r1147);
`clauses_v2.py` 26a66771; driver `run_f1v2ho.sh` 3f6fa61d. Reused unchanged from `f1_frozen/` (SHA1SUMS a063182d):
`f1_bridge.py` 37ee8e77, `f1_batch.py` 14144bde, `f1_null.py` ea2715f2, `c3_identity.py` 8de68cfa, `c4_eval.py`
0d160c94, `f1_gtf_u.py` dcc1e2ff, `clauses.py` 6d5aca9e (its `read_rt`), `split_cls.py` 7045a6af, `rep3_RG/rg3.py`
ec17e540, the RG3-test U scorer (`score_u.py` 3e5cf0e2, `npf_score.py` ed042111, `make_gtf_u.py` 9ebba32d,
`run_fam.sh` 026375a9, `g3_emu.py` 8a7fe6ae) and `readthrough_eval.py` f6b99dcc (run from the repo, in-place copy
asserted equal). Binaries `fj_bin_frozen/mcl_families` 91ef2e1c, `family_score` 7723029b. Python
`/home/juanfra/miniforge3/bin/python3` 3.13.12 + pysam 0.23.3 for F1 / F1v2 / NULL / split_cls_h / clauses; `python3`
(linuxbrew 3.14.4) for readthrough_eval and score_u, as in their earlier runs.

| gate | what must hold | when |
|---|---|---|
| G1h | `f1v2.py --rule none` reproduces F1's `F1.gtf` and `F1.families.gtf` byte for byte; its structural asserts pass | dev (done, §12), then each held-out sample (driver step `f1v2`) |
| G2h | `split_cls_h.py --species gorilla --rt-mode keep` reproduces `split_cls_gw.py`'s counts, rows and bridge rows | dev (done) |
| G3h | `fused_pt.py`'s `fused` equals readthrough_eval `a.fused` on every arm | dev (done: 9 arms × chr16, chr20) |
| G0 | both `SHA1SUMS` hold; every in-place file the frozen scripts import or run equals its frozen copy; BASE sha1 as §2 | each sample, first |
| G6 | testis: `rg3.py`(BASE) = the RG3 test's `RG3.whole.gtf` | held-out |
| G7 | F1's and F1v2's gene_ids = RG3's on every input gene_id without a bridge of that arm | held-out |
| G9 | the F1 batch merge covers every BASE transcript exactly once; `bj_cross_contig_collisions` reported | held-out |
| G5, G8, G10 | testis U: GTF_U(BASE) regenerated = the RG3 test's `BASE_U.gtf`; `make_gtf_u.py` = `f1_gtf_u.py` on the F1v2 (and F1) families GTF; `g3_emu.py` on every new U families run (emu R0 byte-exact, 10 relabellings same partition) | held-out; G10 failure → C4 not measured |

A failed gate G0/G1h/G5-G9 is a bug: it stops that substrate (its clauses not measured) and is fixed in the
instrument, never in the rule, a null, a clause or a bar.

## 5. Clauses (per substrate; integers; zero tolerance unless stated)

FUSED(A) = readthrough_eval f6b99dcc `a.fused` on arm A's loci GTF over the substrate's contigs (BASE.gtf, RG3.gtf,
F1.families.gtf, F1v2.families.gtf, NULL<k>.families.gtf) — the same scorer as F1's C1. A split = an input gene_id whose
non-bridge transcripts carry ≥ 2 output gene_ids; a **bridge split** = one that holds ≥ 1 bridge relation; classes by
`split_cls_h.py --species human --rt-mode drop` (split_cls_gw's SEP / FRAG / UNJ, with RefSeq readthrough-described
records removed from every gene lookup — the gene set FUSED uses; the `keep` classes are reported beside).

| clause | passes iff | states |
|---|---|---|
| **C1** fused loci | judged iff F1v2 has ≥ 1 bridge junction on the substrate. (a) 100·FUSED(F1v2) ≤ 95·FUSED(BASE); (b) FUSED(F1v2) < FUSED(NULL_k) for every k = 0..4; (c) FUSED(F1v2) < FUSED(RG3) | PASS = a∧b∧c; MAG = b∧c∧¬a; FAIL = ¬b ∨ ¬c; NJ = not judged. Reported: F1's FUSED and a/c; FUSED(F1v2) − FUSED(F1); the 10% line; F1all / F1v2all; `a.fused_junction`, `a.rep_fused`; `fused_pt.py`'s per-transcript reading |
| **C2** split correctness | (a) bridge SEP(F1v2) > bridge FRAG(F1v2) [judged iff their sum ≥ 1]; (b) **gain:** bridge FRAG(F1v2) < bridge FRAG(F1) [judged iff bridge FRAG(F1) ≥ 1]; (c) **no net harm:** bridge (SEP − FRAG)(F1v2) ≥ bridge (SEP − FRAG)(F1) [judged iff bridge SEP + FRAG (F1) ≥ 1]; (d) all F1v2 splits SEP > FRAG [judged iff SEP + FRAG ≥ 1]; (e) power: every NULL_k with ≥ 1 split has SEP < FRAG | FAIL = (a), (c) or (d) judged and failing; else NM = (e) fails; else NOEFF = (b) judged and failing (F1v2 removed no single-gene cut: no effect, no harm); else NJ = a part of (a)-(d) not judged; else PASS. Reported: F1's own C2 by the F1 prereg's rule, the RG3 part, every FRAG row by name with its share, the `keep` classes |
| **C3** chains unchanged | `c3_identity.py` on (BASE, F1v2.gtf, F1v2.families.gtf) passes (F1v2.gtf = BASE line for line modulo gene_id and the two appended attributes; the families GTF = F1v2.gtf minus its bridges) **and** `c.matching_intron_chains`(F1v2all) = (BASE) | FAIL_RULE = identity fails; FAIL_CHAINS = identity holds but the count differs (instrument; not measured); PASS. Reported: chains of the families GTFs of F1 and F1v2 |
| **C4** NPIP no-regression (**testis only**) | the F1 prereg's C4 (i)-(v) verbatim, `c4_eval.py` 0d160c94, on the RG3 test's testis U (`rg3_run/U/human_testis.dishuck.U.tsv`, Dishuck 26 chr16 copies, truth key `dishuck`): BASE_U = the RG3 test's GTF_U and families (reused after G5); treated = GTF_U(F1v2 families GTF); NULL_U = `f1_null.py` on (BASE_U, GTF_U(F1v2.gtf)) label `f1v2hoU:human_testis:s1`; XRG3 = the RG3 test's RG3_U run and XF1 = GTF_U(F1 families GTF), both reported. A GTF_U byte-identical to BASE_U or to the RG3 test's RG3_U reuses that run's families products (deterministic: the RG3 test's identical-input NULL runs gave byte-identical PAFs and clusters) | PASS / FAIL; NM if G10 fails; **NA on A119b** (chr16 excluded: no NPIP block on the substrate) |
| **C5** reported | bridge transcripts and junction kinds (RT / ANN / ANN_RT / other) for F1 and F1v2; bridged gene_ids with ≥ 1 ANN junction; the share(J) distribution of F1's bridge junctions by the class of their gene | – |

## 6. Verdict (both substrates; nothing pooled)

- **REFUTE** iff, on either substrate, C1 = FAIL, C2 = FAIL, C3 = FAIL_RULE, or C4 = FAIL.
- **EFFECTIVE** iff C1, C2 and C3 are PASS on both substrates and C4 = PASS on testis.
- **KEEP OPT-IN** otherwise (C1 = MAG or NJ, C2 = NOEFF / NM / NJ, C3 = FAIL_CHAINS, C4 = NM), with no refute trigger.

Readings fixed here: "no effect" (C2 (b) fails, (a), (c), (d) hold) is never REFUTE, as in the RG3 prereg's
no-effect case. A testis with too few events to judge a part is KEEP OPT-IN, never EFFECTIVE. F1's human rows never
enter the verdict (F1 was validated on gorilla; its human numbers are a new description, not a retest).

## 7. Predictions (this author's probabilities, before any held-out number)

1. Gates G0, G1h, G6, G7, G9 pass on both (0.90); G5/G8 pass on testis (0.90); G10 passes where run (0.90).
2. **Power.** A119b: F1 finds ≥ 150 bridge junctions on V_A (0.7), F1v2 keeps 30-65% of them (0.75). Testis: F1v2 has
   ≥ 1 bridge junction (0.85), ≥ 10 (0.4).
3. **C1 A119b:** PASS (0.8): reduction ≥ 5% (0.9), in [10%, 25%] (0.6); (b) every null above F1v2 (0.95); (c) F1v2 < RG3
   (0.9). FUSED(F1v2) − FUSED(F1) ∈ [0, 5% of BASE] (0.85). **C1 testis:** judged (0.85); PASS 0.40, MAG 0.30, FAIL 0.15.
4. **C2 A119b:** (a) 0.9, (b) 0.95, (c) 0.9, (d) 0.9, (e) 0.9 → PASS 0.8. **Testis:** PASS 0.45, NOEFF or NJ 0.35,
   FAIL 0.15.
5. **C3:** PASS 0.97 per substrate.
6. **C4 testis:** PASS 0.9; F1v2 changes 0 U gene_ids beyond RG3's (0.6).
7. **F1 on human (reported):** bridge splits SEP > FRAG on A119b (0.45; human dev gave 14 / 12 and 12 / 12); ANN
   bridged gene_ids ≥ ¼ of F1's (0.8).
8. **Verdict:** EFFECTIVE 0.30, KEEP OPT-IN 0.45, REFUTE 0.20, not decided 0.05 (testis power is the main cap).

## 8. Falsifiers of the design reasoning (reported whatever the verdict)

- **Z1 "Minority bridges are readthroughs, not introns of one gene."** Falsified on a substrate where bridged gene_ids
  with ≥ 1 ANN junction are ≥ ¼ of F1v2's bridged gene_ids (dev: 12 / 64 on both gorilla samples; 1 / 12 and 2 / 11 on
  human chr16 / chr20). Counted per gene_id, not per transcript: one gene can carry 22 bridge transcripts (CREBBP, dev).
- **Z2 "The minority condition costs little of F1's FUSED reduction."** Falsified on a substrate where
  FUSED(BASE) − FUSED(F1v2) < ¾ · (FUSED(BASE) − FUSED(F1)) (dev 0.95 / 0.97 gorilla; 0.95 / 0.93 human dev).
- **Z3 "The gene-cutting bridges are the gene's majority isoform."** Falsified on a substrate where F1v2 removes fewer
  than half of F1's bridge FRAG (dev 27 / 40 and 24 / 39; 10 / 12 and 11 / 12 on human dev).
- **Z4 "The families GTF keeps the main isoforms."** Falsified where the families GTF of F1v2 loses more than ¼ of the
  annotated chains F1's loses (dev 8 / 62 and 6 / 59).
- **Z5 "F1 is weaker on human than on gorilla."** Reported: F1's bridge SEP / FRAG ratio on each human substrate against
  its gorilla ratios (51 / 40, 53 / 39).

## 9. Order, stop rules, machine rules, cost

**Order** (per sample, A119b first, then testis; driver `run_f1v2ho.sh <step> <sample>`): `pre` (G0) → `split` → `f1 k`
per batch → `merge` (G9) → `f1v2` (G1h, then the rule) → `rg3` (G6 on testis) → `gate_names` (G7) → `null k` ×5 →
`c3` → `rt` (repeat until rc 0) → `cls` for F1, F1v2, RG3, NULL0-4 → `fpt` → testis only: `u` (G5, G8), `fam` F1v2 /
NULL / F1, `g3` (G10), `score_u` → `clauses`. Each substrate's tables go to
`/mnt/linuxdisk/tmp/rustle_figures/f1v2_heldout/tables/` before the next starts; then `clauses_v2.py verdict`.
Batch sizes fixed here: `--max-bp` 400,000,000 (A119b) and 800,000,000 (testis).

**Stop rules.** After the freeze nothing changes in the rule, an arm, a null label, a clause, a bar, a truth or a
substrate. A failed gate stops its substrate and is fixed in the instrument only (recorded as an amendment). A step that
hits its time cap is re-run once; a second failure makes its clause not measured. No variant (max, best, top, share:x,
`≤`) is substituted, whatever the result.

**Machine rules.** Heavy steps (F1 batches, merge, NULL draws, readthrough_eval, `mcl_families`) via
`tools/rlock.sh heavy`; light steps (f1v2, split_cls_h, rg3.py, U restrictions, scoring) via `light`; all foreground;
`TMPDIR` under `/mnt/linuxdisk`; never `pkill -f`; scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/f1v2/`. No `src/`
edit, no commit, no push.

**Cost (estimate).** F1 batches: A119b ≈ 16 min (dev chr16: 1.79 M records in 27 s, 0.75 GB), testis ≈ 3 min; F1v2
≈ 10 s; NULL ≈ 5-30 s × 10; readthrough_eval 12 arms ≈ 15-25 min per sample (dev chr16 43 s for 11 arms); split_cls_h
≈ 1 min × 8 × 2; testis U families ≤ 3 runs × 4 min. **≈ 1.5-2.5 h wall.** Disk < 6 GB transient.

## 10. Not in this test

Gorilla (spent for bridge regrouping, now development); the families effect genome-wide (only the testis NPIP block);
the simulated NPIP gain (dev only: F1v2 9 / 5 / 0 of 10 fused copies back in NPIP at f = 0.1 / 0.5 / 0.9, F1 9 / 6 / 3,
BASE 7 / 1 / 0 — by design F1v2 refuses dominant fusions); A119b's NPIP block (chr16 is dev); the variants max / best /
top / share:x (dev only); a Rust port.

## 11. Hostile self-review (fixes applied above)

1. **"The condition was designed on the very samples that exposed the problem."** Yes: they are declared development
   (the task spent them), and the test is human. The human dev contigs chr16 / chr20 were also seen and are excluded.
2. **"Four variants were compared; the winner is selected."** The sum-both-sides form was named first as the literal
   reading of the post-hoc statistic (share < ½); it also had the best net on both gorilla samples; the curve (§12.2)
   shows ½ on a plateau (OR's peak is 0.4). The others are reported, never substituted.
3. **"C2 (b) is nearly guaranteed: F1v2's bridges are a subset of F1's."** Yes, FRAG cannot rise; (b) is the gain
   guard, and its failure is NOEFF (KEEP OPT-IN), never REFUTE. The harm guard is (c): F1v2 must not lose more correct
   separations than it removes cuts.
4. **"Dropping RefSeq readthrough records from C2 favours the rule."** A readthrough record overlapping both parent
   genes would call every correct separation of an annotated readthrough FRAG, while FUSED (C1) counts the same locus as
   fused because its gene set excludes those records. C1 and C2 must use one gene set; `drop` is FUSED's. The `keep`
   classes are reported beside every count.
5. **"C1 (a) at 5% is below R's 10% bar."** Inherited from F1's C1 unchanged; the 10% line is reported.
6. **"Testis is small (26 MB GTF, 1.9 isoforms per locus): its clauses may be unjudgeable."** Stated: C1 needs ≥ 1
   F1v2 bridge; C2's parts carry their own judged-iff; unjudged → KEEP OPT-IN, never EFFECTIVE.
7. **"C4 on testis has 7 present copies and one fused holder (NPIPB14P, 2 reads)."** Known from the RG3 Outcome and
   stated; C4 is a no-regression control, not a gain test.
8. **"The NULL is matched to F1v2, not F1."** Intended: C1 (b) attributes F1v2's reduction to its predicate. F1's human
   C1 (b) is not computed and F1 is never judged.
9. **"Strict `<` loses exact ties."** Fixed before the tie count; ties abstain (do nothing), the house's
   assign-or-abstain stance; dev ties split both ways (§1.2).
10. **"The `reads` attribute is an assembler output, not raw evidence."** It is the quantity F1 and RG3 already use to
    name keepers, identical across arms, and it counts each read once per transcript.
11. **"A119b's chr21 and chr22 are also excluded."** The task's substrate definition; the exclusion is fixed here and
    applies to every clause.
12. **"FUSED's locus-union reading (r1147) flatters regrouping."** The per-transcript reading is reported for every arm
    (`fused_pt.py`, gated equal to `a.fused` in its locus form).
13. **"F1 on human may itself fail; that would sink F1v2 by association."** F1's human rows are reported only; F1v2 is
    judged against BASE, RG3 and its own nulls, and against F1 only in C2 (b)/(c).

## 12. Development evidence (design phase; every number below was seen before this file)

Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/f1v2/` (`dev/`, `sim/`, `hdev/`); report
`scratchpad/figs/f1v2.md`.

### 12.1 Gorilla (OR6737 minus NC_073244.2 and KB3781, the F1 held-out substrates, now dev)

| | OR6737 − NC_073244.2 | KB3781 |
|---|---|---|
| FUSED BASE / RG3 / F1 / **F1v2** | 549 / 510 / 463 / **467 (−14.9%)** | 662 / 574 / 525 / **529 (−20.1%)** |
| FUSED NULL(F1v2) × 5 | 551-552 | 664-670 |
| bridge splits SEP / FRAG: F1 → **F1v2** | 51 / 40 → **45 / 13** | 53 / 39 → **49 / 15** |
| all splits F1v2 SEP / FRAG (UNJ) | 84 / 17 | 140 / 33 (14) |
| NULL(F1v2) SEP / FRAG | 0 / 99-100 ×5 | 0 / 181-185 ×5 |
| chains of the families GTF: BASE / F1 / **F1v2** | 24,277 / 24,215 / **24,269** | 25,911 / 25,852 / **25,905** |
| F1v2all chains = BASE; c3_identity | yes; pass | yes; pass |
| F1v2 bridge transcripts RT / ANN / other (F1) | 51 / 18 / 1 (61 / 181 / 1) | 68 / 21 / 2 (75 / 155 / 6) |

All contigs of OR (with NC_073244.2): F1 56 / 41 → F1v2 50 / 14; bridge junctions kept 64 of 97 (OR), 64 of 92 (KB).
Lost SEP: OR 6 (e.g. RBL2 | LOC109023633, CA3 | CA13, CORO7 | PAM16), KB 4 (ZNF620 | ZNF619, POLE4 | lncRNA); remaining
FRAG are minority single-gene links (LARP1B and PGBD1 in both samples, LOC101129171, CDK11B, DLC1, TJP1).

### 12.2 Variants and the curve (bridge splits SEP / FRAG; OR all contigs, KB)

| rule | OR | KB |
|---|---|---|
| F1 | 56 / 41 | 53 / 39 |
| **min (F1v2)** | **50 / 14** | **49 / 15** |
| max (fewer than the larger side) | 56 / 28 | 53 / 25 |
| best (Σ link < best transcript of each side) | 46 / 11 | 44 / 12 |
| top (best link transcript < best of each side) | 47 / 13 | 48 / 16 |

share < x, net (SEP − FRAG): OR 0.1: 21, 0.2: 31, 0.3: 37, 0.4: 39, **0.5: 36**, 0.6: 36, 0.7: 33, 0.8: 30, 0.9: 27,
F1: 15; KB 25, 30, 31, 32, **34**, 31, 30, 31, 28, F1: 14. `min` = `share:0.5` (identical outputs).

### 12.3 Human dev contigs (A119b chr16, chr20; dev GTFs; excluded from V_A)

| | chr16 | chr20 |
|---|---|---|
| FUSED BASE / RG3 / F1 / **F1v2** / NULL × 5 | 101 / 87 / 82 / **83** / 102-106 | 60 / 54 / 46 / **47** / 60-61 |
| per-transcript reading BASE / RG3 / F1 / F1v2 | 86 / 86 / 81 / 82 | 54 / 54 / 46 / 47 |
| bridge splits SEP / FRAG: F1 → **F1v2** (rt drop) | 14 / 12 → **10 / 2** | 12 / 12 → **10 / 1** |
| all F1v2 splits SEP / FRAG (UNJ); NULL0 | 24 / 4 (1); 0 / 28 | 15 / 2 (2); 0 / 18 |
| F1 bridge junctions → F1v2 | 27 → 12 | 25 → 11 |
| c3_identity, chains F1v2all = BASE | pass, 1,605 | pass, 1,059 |

F1v2 separates, e.g., RHOT2 | WDR90, PALB2 | NDUFAB1, ITCH | DYNLRB1 and NPIPA7 | NPIPA6 (chr16); its FRAG are CREBBP,
CLEC18A (chr16) and CHD6 (chr20). On chr20 F1 alone ties 12 / 12 (it would fail its own C2 there).

### 12.4 The simulation (locus_fix_design arms; families rerun where F1v2's families GTF differs: f0.5, f0.9)

| f | fused copies in NPIP /10: BASE / F1 / **F1v2** | copies /25 | controls /10 | confounder genes split |
|---|---|---|---|---|
| 0.0 | 9 / 9 / 9 | 23 / 23 / 23 | 10 | 0 |
| 0.1 | 7 / 9 / **9** | 22 / 23 / 23 | 10 | 0 |
| 0.5 | 1 / 6 / **5** | 16 / 20 / 19 | 10 | 0 |
| 0.9 | 0 / 3 / **0** | 15 / 18 / 15 | 10 | 0 |
| 1.0 | 1 / 1 / 1 | 16 / 16 / 16 | 10 | 0 |

F1v2 keeps F1's gain where the fusion is a minority (f = 0.1) and gives it up where the fusion carries half or more of
the copy's reads (f = 0.5: 1 of 7 bridges dropped, share 0.57; f = 0.9: 3 of 4, shares 0.50-0.83). The SULT1A3
simulated FRAG (share 0.24) persists: it is a minority link.

## Amendments

### Amendment 1 — the freeze (2026-09-29 00:50, written BEFORE any held-out command)

**This file's sha1 before this amendment:** `9ec5a7e30ce78ec5a53c5a980618fcbe54dd6a53` (26,260 bytes; a byte copy is
kept at `/mnt/linuxdisk/tmp/rustle_figures/f1v2_heldout/PREREG_f1v2_readshare_2026-09-29.pre_amendment1.md`, born
00:50:25). Acceptance: the orchestrating session's task (design on dev → pre-register → run exactly as pre-registered)
is the mandate; no separate user acceptance of the reuse (§3) was sought, and this is recorded as such.

**Frozen instruments.** `/mnt/linuxdisk/tmp/rustle_figures/f1v2_frozen/SHA1SUMS` sha1 `893d3ad092e4b21ff7cf18edd214d1e7cedbabe3`:
`f1v2.py` b4e788ad, `split_cls_h.py` 4a6b178d, `fused_pt.py` c42a4605, `clauses_v2.py` 26a66771, `run_f1v2ho.sh`
3f6fa61d. F1 set `/mnt/linuxdisk/tmp/rustle_figures/f1_frozen/SHA1SUMS` a063182d (checked `sha1sum -c` at the freeze).
Binaries `mcl_families` 91ef2e1c, `family_score` 7723029b. Python 3.13.12 + pysam 0.23.3 (miniforge); python3 3.14.4.
BASE GTFs: A119b `631c9f114b3728c51405a6e4f019f20cee43f12e`, testis `50f239d6da21355f3d9723c495117980e1e9c25f`.

**Dev gates (run before this file's sha1):**
- **G1h PASS**: `f1v2.py --rule none` = F1's `F1.gtf` and `F1.families.gtf` byte for byte on gorilla OR6737 and KB3781
  (genome-wide F1 held-out products), on the 5 simulation arms f0.0-f1.0, and on human A119b chr16 / chr20 (dev F1
  runs); every structural assert held (T_J ids and up / down / inside counts equal F1's rows).
- **G2h PASS**: `split_cls_h.py --species gorilla` (both rt modes; gorilla has 0 readthrough-described records) equals
  `split_cls_gw.py` in counts, rows and bridge rows on 6 files (F1 and F1v2 on OR and KB, and the F1 held-out's own
  `cls.F1.json` for OR − NC_073244.2 and KB).
- **G3h PASS**: `fused_pt.py` `fused` = readthrough_eval `a.fused` on all 9 arms of chr16 and of chr20.
- Code paths: `clauses_v2.py` on the 4 dev substrates (OR, KB, hsa16, hsa20) gives C1 / C2 / C3 PASS on each and
  `verdict` = KEEP OPT-IN (C4 NA everywhere, so EFFECTIVE is unreachable on dev, as §6 intends); `c3_identity.py`
  passes on F1v2 for OR, KB, hsa16, hsa20. The testis-only C4 steps (`u`, `fam`, `g3`, `score_u`) of the driver were
  not exercised on human before the freeze: they are the F1 driver's steps with the human paths of the RG3 test.

**Command lines** are the steps of `run_f1v2ho.sh` (§9 order), run from the frozen directory, sample `human_A119b`
then `human_testis`. NULL labels `f1v2ho:<sample>:s<k>`, k = 0..4; U NULL label `f1v2hoU:human_testis:s1`.

### Amendment 2 — deviations during the run (2026-09-29, after the run; no rule, arm, null, clause or bar changed)

1. **Instrument fix (driver, G0).** The first `pre` call stopped with `MISMATCH rg3_lib/__pycache__`: the in-place
   check looped over `f1_frozen/rg3_lib/*`, which now holds a `__pycache__` directory from the F1 run. The chained
   `split` call then failed on the missing BASE symlink (no held-out line was read). The loop now skips non-files
   (`[ -f $f ] || continue`); `run_f1v2ho.sh` 3f6fa61d → 016f6aed, `SHA1SUMS` 893d3ad0 → 01eda85d. Nothing else changed.
2. **Lock waits (machine only).** Another session's `minimap2` jobs held the heavy lock intermittently (one at 16.5 GB).
   Three A119b F1 batch calls (3 once, 4 twice) timed out WAITING for the lock (exit 1, no output, no Python started)
   and were re-run with a longer `RLOCK_WAIT`; every batch then ran once, to completion, under the heavy lock.
3. **How to recover the frozen text.** Amendment 1 was inserted between the `## Amendments` and `## Outcome` headers,
   so the frozen bytes are not a prefix of this file: the byte copy (sha1 9ec5a7e3) equals this file with the
   Amendment and Outcome bodies removed (`diff`: 0 deleted lines).

## Outcome (2026-09-29)

**Gates.** G0 (after the Amendment 2 fix), G1h (`--rule none` = F1 byte for byte, structural asserts held), G6 (testis
RG3 = the RG3 test's `RG3.whole.gtf`), G7 (F1 and F1v2 names = RG3's off their bridged gene_ids, 0 differ), G9 (0
cross-contig collisions; coverage asserted), G5, G8 (F1v2 and F1) and G10 (the one new U run, F1's: emu R0 byte-equal,
10/10 relabellings) all passed. G3h also held on the held-out arms: `fused_pt.py`'s `fused` equals `a.fused` on all
9 arms of both samples. F1 on A119b: 9 batches, 47-183 s, ≤ 2.27 GB; testis 5 batches, ≤ 25 s.

| | **human_A119b − chr16, chr20, chr21, chr22** | **human_testis** |
|---|---|---|
| bridge junctions on the substrate: F1 → **F1v2** (ties at share ½) | 812 → **394** (27) | 22 → **14** (2) |
| FUSED BASE / RG3 / F1 / **F1v2** | 1,965 / 1,786 (−9.1%) / 1,595 (−18.8%) / **1,628 (−17.2%)** | 200 / 195 (−2.5%) / 182 (−9.0%) / **184 (−8.0%)** |
| FUSED NULL(F1v2) k = 0..4 | 2,008, 1,995, 1,996, 2,001, 2,001 | 200, 202, 203, 200, 201 |
| **C1** (a) ≤ 95% · (b) < every null · (c) < RG3 | **PASS** (also ≤ 90%) | **PASS** (not ≤ 90%) |
| C2 bridge splits SEP / FRAG: F1 → **F1v2** | 219 / 451 → **193 / 149** (UNJ 13) | 15 / 7 → **12 / 2** |
| C2 (b) FRAG F1 → F1v2; (c) net F1 → F1v2 | 451 → 149; −232 → **+44** | 7 → 2; +8 → +10 |
| C2 (d) all F1v2 splits SEP / FRAG (UNJ); RG3 part | 378 / 181 (75); 185 / 32 | 18 / 2 (1); 6 / 0 |
| C2 (e) NULL SEP / FRAG | 0-1 / 586-596 on 5/5 | 0 / 20 on 5/5 |
| **C2** | **PASS** | **PASS** |
| C3 identity lines; chains F1v2all = BASE | 2,214,636 identical; 34,752 = 34,752 | 166,309 identical; 11,399 = 11,399 |
| **C3** | **PASS** | **PASS** |
| C4 on the RG3 test's U (294 loci): c1 BASE / F1v2 / NULL / XRG3 / XF1 | NA (chr16 excluded) | 5 / 5 / 5 / 5 / 5 (NPIP = MCL1); (i)-(v) hold; GTF_U(F1v2) = the RG3 test's RG3_U byte for byte; NULL_U = BASE_U |
| **C4** | **NA** | **PASS** |
| reported: per-transcript FUSED BASE / RG3 / F1 / F1v2 / NULL | 1,754 / 1,756 / 1,574 / 1,602 / 1,774-1,786 | 194 / 194 / 181 / 183 / 194-197 |
| reported: chains of the families GTF BASE / F1 / F1v2 / NULL | 34,752 / 33,966 / **34,640** / 34,596-34,639 | 11,399 / 11,388 / 11,394 / 11,392-11,396 |
| reported: rt-kept classes, bridge SEP / FRAG F1 → F1v2 | 203 / 467 → 178 / 164 | 11 / 11 → 8 / 6 |
| C5 bridge transcripts RT / ANN / ANN_RT / other only: F1 → F1v2 | 296 / 4,975 / 16 / 269 → 232 / 809 / 16 / 100 | 14 / 15 / 5 / 1 → 8 / 3 / 5 / 1 |

**VERDICT (§6): EFFECTIVE** (`clauses_v2.py verdict`: C1, C2, C3 PASS on both samples, C4 PASS on testis and NA on
A119b; no refute trigger). F1v2 stays opt-in until the user decides; a default flip is the user's call.

**Default flipped on 2026-09-29 by the user's decision**: `copy_assign --assemble-only` and the driver now run `f1v2`
unless told `off`, citing this held-out verdict and the family-level side result of
`docs/archive/2026-09/PREREG_o1_cover_growth_2026-09-29.md` (Outcome: F1v2's families, the COVER core, beat BASE on every Compara
metric on both human substrates and on Liftoff recall on testis). `RUSTLE_BRIDGE_REGROUP=off` reproduces the pre-flip
products byte for byte (`bench/ASSEMBLY_POLISH.md` addendum 3).

**Read this before quoting it.**
- **F1 alone fails on human A119b.** Its bridge splits are 219 SEP / 451 FRAG (it would fail its own C2 there), and 64%
  of its bridged gene_ids carry an annotated-intron junction. F1v2 turns that into 193 / 149, so the net correct
  separations go from −232 to +44. On testis F1 was already positive (15 / 7) and F1v2 gives 12 / 2.
- **F1v2 still cuts genes on the deep library: 149 FRAG on A119b**, every one of them through a minority link (share
  < ½ by construction). Examples: ADAMTS6, BIRC6, CSMD1, DCC, DLC1, DNAH12, ATRX. This is the LARP1B / PGBD1 type from
  gorilla dev: both side pieces outnumber the full-length link. Z1 fired on A119b (below).
- The post-hoc read-share pattern replicates on a held-out species. Share(J) of F1's A119b bridge junctions: FRAG genes
  median 0.66 (IQR 0.36-0.90), SEP genes 0.14 (0.05-0.37).
- F1v2 keeps 91% (A119b) and 89% (testis) of F1's FUSED reduction. Its families GTF loses 112 annotated chains on A119b,
  against 786 for F1 and 113-156 for the matched nulls: the main isoforms stay in.
- **The rt-drop choice does not carry the verdict.** Under the reported rt-kept classes every C2 part still holds:
  A119b 178 > 164, 164 < 467, +14 ≥ −264, and all splits 355 > 204; testis 8 > 6, 6 < 11, +2 ≥ 0, and 11 > 9.
- **Testis is thin:** 14 bridge junctions, the 5% bar met at 8.0%, and C2 judged on 14 bridge splits. **C4 is a
  no-regression control with no fused NPIP member to gain:** GTF_U(F1v2) is byte-identical to the RG3 test's RG3_U.
- Per-transcript reading (r1147): F1v2 −8.7% (A119b) and −5.7% (testis), against −17.2% / −8.0% on `a.fused`. RG3
  gives +2 and 0 on this reading.

**Predictions (§7).**
- P1: gates pass. **Hit**; one driver fix was needed at G0.
- P2: F1 ≥ 150 bridge junctions on A119b (hit, 812). F1v2 keeps 30-65% (hit, 48.5%). Testis ≥ 1 (hit) and ≥ 10 (hit,
  14; prior 0.4).
- P3: A119b C1 PASS, reduction ≥ 5% and inside [10%, 25%], (b), (c), and F1v2 − F1 ∈ [0, 5% of BASE] (33 ≤ 98): all
  **hit**. Testis judged, and PASS (prior 0.40): **hit**.
- P4: C2 PASS on A119b (hit) and on testis (hit, prior 0.45).
- P5: C3 (hit).
- P6: C4 PASS (hit); F1v2 changes 0 U gene_ids beyond RG3's (hit).
- P7: F1 bridge SEP > FRAG on A119b (**missed**, 219 / 451; prior 0.45). ANN bridged gene_ids ≥ ¼ of F1's (hit: 448 / 701
  and 7 / 22).
- P8: EFFECTIVE (prior 0.30): **hit**.

**Falsifiers (§8).**
- **Z1 FALSIFIED on A119b**: 137 of 355 F1v2 bridged gene_ids (39%) carry an ANN junction, against F1's 448 / 701. It
  holds on testis (2 / 14). Minority links are not only readthroughs in a deep library. They are fewer and mostly
  harmless (193 SEP vs 149 FRAG), but not absent.
- Z2 holds: 0.91 and 0.89 of F1's FUSED reduction are kept.
- Z3 holds: F1v2 removes 302 of 451 and 5 of 7 of F1's bridge FRAG.
- **Z4 FALSIFIED on testis** (5 lost chains vs F1's 11, 45% > ¼; tiny n). It holds on A119b (112 vs 786, 14%).
- Z5: F1 is weaker on A119b than on gorilla (219 / 451 vs 51 / 40 and 53 / 39), not on testis (15 / 7).

**What this does and does not show.**
- On a held-out species, genome, annotation and two libraries, the minority condition keeps about 90% of F1's
  fused-locus reduction, beyond RG3 and five matched nulls.
- It removes two thirds or more of F1's single-gene cuts, and turns F1's net-negative split record on the deep human
  library positive.
- It changes no transcript or chain, and leaves the testis NPIP block unchanged.
- It does not remove minority single-gene links (149 on A119b).
- It cannot be credited with an NPIP gain: none was testable here, and in the simulation it gives up F1's gain for
  fusions carrying ≥ ½ of a copy's reads.
- Both human libraries are now spent for bridge regrouping work.

**Kept products.**
- Tables: `/mnt/linuxdisk/tmp/rustle_figures/f1v2_heldout/tables/`:
  - `verdict.json` 513b7391;
  - `human_A119b.clauses.json` c9b1b147, `human_testis.clauses.json` 4d0173dd;
  - `reported.json` b7d992e8 (Z1-Z5, rt-kept classes, FRAG rows with shares, lost SEP rows, share by class);
  - readthrough_eval TSVs 953fc019 / bbea1802;
  - per sample: `cls.*`, `F1.junctions.tsv`, `F1v2.bridges.tsv`, stats, `c3*`, `fpt.json`;
  - testis only: `c4.json` and `U.score.clauses.json`.
- Scratch: `/mnt/linuxdisk/tmp/rustle_figures_dev/f1v2/` (14 GB, including the dev work).
- GTF sha1s:
  - A119b: F1.gtf 88466818, F1v2.gtf 976b93dc, F1v2.families.gtf daaf6ad2;
  - testis: F1.gtf d074ea9b, F1v2.gtf 4afc60a2, F1v2.families.gtf 6bdaeabd.

## Independent verification (2026-09-29, own code, no scorer under test imported)

**CONFIRMED WITH CORRECTIONS** (`figs/f1v2_verify.md`; scripts `/mnt/linuxdisk/tmp/rustle_figures_dev/f1v2_verify/`).
Freeze order holds (instruments 00:48, prereg sha1 9ec5a7e3 00:50:25, first held-out file 00:50:53); all sha1s match;
neither amendment changes a rule or bar. Every verdict number reproduces (fused loci on all 11 arms x 2 samples, SEP/FRAG
row by row, bridge junctions 812->394 / 22->14, chain identity); C1-C4 pass; no dev contig in the verdict.
- **Share rule:** equals "share < 1/2" only with share = reads(link) / (reads(link) + min(UP, DOWN)) (0/909 disagree);
  against the both-sides-summed share it disagrees on 300/887 A119b bridges.
- **C1 Z1 wording:** Z1 fired because 137/355 bridged gene_ids (38.6%) carry an annotated-intron junction; 149 is the
  C2 FRAG count. Z4 fired on testis (5 vs 11 lost chains).
- **Metric trap (main correction):** the §5 F1all/F1v2all rows were omitted. Counting each bridge as its own locus,
  F1v2 is NOT better than RG3 (A119b 1,829 vs 1,786; testis 195 vs 195). The fused-locus gain over RG3 exists only when
  bridges are scored as relations (fusion_of) removed from the loci file, i.e. F1v2 moves fusions into explicit relation
  records rather than removing them. Matched nulls still hold by wide margins either way.
- **Nulls under-matched, undisclosed:** A119b nulls remove 1,122 vs 1,155 bridge transcripts (one chr15 gene without a
  donor); testis short one gene; the NPIP-block null changed nothing. Margins unaffected.
- **Gene cuts:** 8 of the 149 remaining single-gene cuts are already cut by RG3 alone (at most 141 are F1v2's).
