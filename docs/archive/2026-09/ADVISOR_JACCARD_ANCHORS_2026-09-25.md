# Advisor's proposal: read-Jaccard + multimapper chaining anchors instead of exon sums

## §0 Pre-registration (written 2026-09-25 before any arm was run)

**Proposal (Canzar, via the user):** comparing exon sums between loci is overcomplicated. Instead, use the
**Jaccard index of the read sets** of two loci, and use the multi-mapping reads both for **connectivity**
(a shared molecule links two paralogous loci) and as **chaining anchors** (the same read base aligned at
locus A and at locus B is an anchor (posA, posB); chain the anchors collinearly and measure how much of each
locus the chain covers).

**Instrument** (`bench/mechanism/jaccard_anchor_test.py`). GUIDED nodes: every annotated gene/pseudogene on
the contig with its exon union, so every truth gene is a node in every arm and **only the edge rule changes**.

| arm | edge rule |
|---|---|
| `S_shipped` | shipped: `minimap2 -x asm20` all-vs-all of gene spans → `mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60` (its own MCL) |
| `S_port` | the same pre-MCL graph (`--dump-graph`) through `mcl_port` I=2.8 — the operator control for the read arms |
| `J(t)` | Jaccard of molecule sets (primary **or** secondary block with ≥20 bp on the gene's exons) ≥ t, t ∈ {.01,.05,.10,.20,.30}, weight = J |
| `C` | multimapper anchors, LIS-chained, edge iff the chain covers ≥ 0.60 of the SMALLER gene's exonic bases (the exon-sum rule rebuilt from reads, same C) |
| `J0.05+C` | both (the advisor's full combination) |
| `S\|C` | union — does the read chain add pairs the DNA rule misses? |

Gene pairs whose exon unions overlap (one locus) are never read edges. Read arms cluster with `mcl_port` I=2.8.
Scored by `family_score --pairs` (sensitivity / precision / bipartite F / collapsed + pairwise) against the
protein referee, plus edge-level TP/FP and a **ceiling table** (what fraction of truth same-family pairs share
≥1 molecule at all).

**Substrates.** DEV: human A119b chr20 (1.10 M records, 680 k secondary; referee 42 fams / 152 genes).
HELD-OUT: gorilla testis `GGO_mm` NC_073244.2 (473 k records, 314 k secondary; referee 113 fams / 775 genes).
Both BAMs are `-N 50 -p 0.1 --secondary=yes`. Never pooled. The Jaccard threshold t is picked on DEV (best F)
and applied unchanged to held-out.

**Decision rule (held-out gorilla, bipartite F vs `S_shipped`):**
- ⭐ **WORKS** — `J(t*)`, `C` or `J0.05+C` has F ≥ S_shipped + 0.02 AND precision not lower by > 0.05.
- ⚠ **EQUIVALENT** — |ΔF| < 0.02: the simpler rule could replace exon sums; report it as a simplification.
- ⛔ **DOES NOT WORK** — F lower by ≥ 0.02.
- `S|C` is reported separately: if it beats S_shipped by ≥ 0.02, reads are a **complement**, not a replacement.

**Predicted before looking:** ⛔ for J and C as replacements. (1) A read rule can only see pairs where BOTH genes
are expressed in the library and minimap2 emits the secondary — r1026 already bounds even genomic alignment at
~19% of referee pairs, and read secondaries need far higher identity than asm20. (2) Jaccard is dominated by
expression imbalance (a 1000-read parent vs a 5-read paralog has J ≤ 0.005 even if every paralog read is
shared) — register 293/359 already refuted minimizer Jaccard for the length-mismatch analogue. (3) Shared
molecules also link non-homologous genes through repeats in UTRs (register 721: 27 of 92 read-linked family
pairs had no homology). Expect `S|C` ≈ S (the read-linked pairs are the high-identity ones the DNA rule already
has; r1059's ALL-arm gain was connectivity from the same chains).

### §0b Addendum (written after dev + gorilla were seen, BEFORE chr16 was run)

Two things were added after the dev run and are **not** covered by §0:
1. **Multimapper-only sharing.** On dev, 14 of 30 read-linked truth pairs were linked only by ONE alignment
   spanning both genes (readthrough between tandem neighbours), which is not the advisor's multimapping idea. So
   `J(t)` now counts a molecule only if two DIFFERENT records reach g and h, neither spanning both. The old
   count is kept as `Jany(t)` (a control).
2. **`Cany` — the anchor chain as a CERTIFICATE, not a coverage floor.** Keep a multimapper link only if the two
   alignments share ≥1 read base (cov > 0). Motivation: 249/389 dev multimapper pairs had cov = 0 (the read's
   two alignments cover DISJOINT read segments, so they are split alignments, not paralogous placements).

Because `Cany` was chosen after seeing the numbers, it is tested on a **third substrate (human A119b chr16,
NPIP/LCR16 rich, referee `/mnt/linuxdisk/tmp/sedef/truth/chr16.tsv`)** with the arms frozen: `J0.01`,
`Cany`, `S_shipped`, all `.cc`. Same §0 bar: F ≥ S_shipped + 0.02 AND precision within 0.05.

## §1 Verdict

**Multimapper Jaccard is NOT a drop-in replacement for the exon-sum rule, and the chaining-anchor version of the
exon-sum rule is strictly worse.** Reads link many more real paralog pairs than the DNA rule (2–7× the truth
pairs), but they also link genes that were **duplicated together without being homologous** (segmental-
duplication passengers, chimeric models, a gene that carries a piece of another's exon). The exon-sum fraction
is exactly the test that rejects those links. On the young-SD substrate (chr16), family precision drops from .96
to .72. Read-built chaining with the same exon-sum threshold keeps the precision but loses the recall, because
reads never cover a gene's whole exon union.

| arm (held-out gorilla, pre-registered) | F | prec | ΔF vs `S_shipped` | Δprec | verdict (§0 bar) |
|---|---|---|---|---|---|
| `S_shipped` | .208 | .989 | — | — | — |
| `J0.01` (t* from dev) | .455 | .849 | **+.247** | **−.140** | ⛔ fails the precision clause: a recall-for-precision trade, not "works" |
| `C` (anchor chain ≥ .60 of smaller gene) | .178 | .987 | −.030 | −.002 | ⛔ **does not work** |
| `J0.05+C` (the full combination) | .178 | .987 | −.030 | −.002 | ⛔ **does not work** |
| `S\|C` (union) | .281 | .992 | +.073 | +.003 | ⭐ passes as a **complement**, but see the operator control below (+.033) |

## §2 All three substrates (`.cc` = connected components, the operator-free control; `S.cc` is the fair baseline)

Family level, bipartite (`family_score`): sensitivity / precision / **F**.

| arm | DEV human chr20 (42 fams) | HELD-OUT gorilla NC_073244.2 (113 fams) | 3rd: human chr16 (LCR16/NPIP) |
|---|---|---|---|
| `S_shipped` (its own MCL) | .092 / 1.000 / **.169** | .116 / .989 / **.208** | .176 / .964 / **.298** |
| `S.cc` (same graph, components) | .132 / .952 / **.231** | .142 / .991 / **.248** | .232 / .973 / **.375** |
| `J0.01` multimapper Jaccard | .164 / .862 / **.276** | .311 / .849 / **.455** | .258 / .725 / **.381** |
| `J0.05` | .151 / .920 / .260 | .246 / .960 / .392 | .252 / .819 / .385 |
| `J0.30` | .112 / .895 / .199 | .186 / .966 / .312 | .248 / .884 / .388 |
| `Jany0.01` (incl. readthrough) | .342 / .852 / .488 | .364 / .803 / .501 | .363 / .730 / .485 |
| `C` chain ≥ .60 | .013 / 1.000 / .026 | .098 / .987 / .178 | .069 / .955 / .128 |
| `C0.30` chain ≥ .30 (descriptive) | .026 / 1.000 / .051 | .128 / .980 / .226 | .176 / .964 / .298 |
| `Cany` chain > 0 (post hoc) | .118 / .900 / .209 | .307 / .971 / **.467** | .219 / .713 / .335 |
| `S\|C` | .132 / .952 / .231 | .164 / .992 / .281 | .239 / .973 / .383 |

`.mcl` arms (`mcl_port` I=2.8) are in `<tag>.summary.tsv`. At I=2.8 MCL breaks the sparse read graphs into
singletons, so J scores the same at every t on dev. Components are the literal reading of "threshold the Jaccard
and connect".

**Ceiling: which truth same-family pairs can a read rule see at all?**

| | chr20 | gorilla | chr16 |
|---|---|---|---|
| truth pairs | 404 | 41,515 | 1,474 |
| both genes expressed (≥3 primary molecules) | 72.8% | 82.2% | 88.9% |
| ≥1 **multimapping** molecule | **4.0%** | **8.1%** | **20.5%** |
| shipped DNA edge | 2.0% | 1.2% | 10.6% |
| multimapper link but no DNA edge | 9 | 2,905 | 149 |
| multimapper links between DIFFERENT truth families | 7 | 237 | **509** |

## §3 Why: the four mechanisms

1. **Reads find paralogs the DNA rule misses, especially retrocopies.** The gorilla families that go from F 0 to
   F 1 are almost all a parent with its processed copy (FTL~LOC101133655, RPS15~LOC101129785, RPL36~LOC109024144,
   APOC1, TPRX1/2). A spliced read aligns contiguously to an intronless retrocopy, so **a read already is an exon
   sum**. The DNA rule has to rebuild that sum from a genomic alignment that introns break up. That is the true
   part of the advisor's intuition.
2. **Reads also link co-duplicated passengers.** A segmental duplication copies a block holding several genes
   and pieces of genes. A read from one copy multimaps to the other copy of the block and lands on whatever
   gene models sit there. chr16's top false links: **EIF3CL~NPIPB9 (5,424 molecules, J = .47)**, where EIF3C
   carries ~0.8 kb of the LCR16a/NPIP core (known, `mcl_families` `--core-refine` doc); and the BOLA2 /
   BOLA2-SMG1P6 / SLX1A / LOC… block (the 16p11.2 BP4–BP5 duplicon). Jaccard only measures **"these loci were
   duplicated together and are expressed"**. On a young SD it cannot tell a family from a duplicon. This is
   register rows 721 (27 of 92 read-linked family pairs had no homology), 1081 and 1190, now measured as a
   family-level precision loss. **The exon-sum fraction is the component that rejects these links**, which is
   why it is there. It is not an over-complication.
3. **Jaccard is expression-weighted, so t is not a similarity.** One molecule set is over the union. A
   1,000-read parent and a 5-read paralog top out at J ≈ .005 even if every paralog read is shared. So the best
   t is tiny (.01, dev and held-out). Raising t removes true pairs from unbalanced families faster than it
   removes co-duplication links (chr16 precision only reaches .88 at t = .30, still below the DNA rule's .96).
   Canzar's own objection to arbitrary similarity thresholds applies here, and t moves with library depth and
   tissue.
4. **Chaining the anchors cannot rebuild the exon sum from reads.** The anchors are real (a shared read base
   aligned at both loci), but reads cover only the expressed part of each gene, and `-p 0.1` secondaries are
   often local. The chain covers ≥ .60 of the smaller gene's exon union in only 210/3,343 true gorilla pairs
   (6%). At that threshold `C` is precise but sees almost nothing. Relaxed to .30 it only **matches** the
   shipped rule on chr16 (.298 = .298) and loses on the other two. Separately, 41% of true and 79% of false
   gorilla multimapper links have **zero** shared read bases: the read's two alignments are disjoint pieces of
   the read (split/partial secondaries), so they are not paralogous placements. That is the post hoc `Cany`
   certificate. It helped on gorilla (F .467, prec .971) but **failed its frozen third-substrate test** (chr16
   precision .713), because co-duplicated passengers DO share read bases.

## §4 What this means for the thesis

- **Keep the exon-sum rule as the definition.** Its job is to decide whether two loci share enough of their
  exons, as opposed to sharing a duplicated segment. Reads cannot make that decision: they see the duplicated
  segment only when it is expressed, and they cannot see the part of the exon union that is not transcribed.
- **Architectural cost (unchanged from r1081):** defining families by multimapping makes O1 depend on the
  aligner's ambiguity, which is O2's input. Jaccard is a support relation, not an assignment, so the O1⊥O2
  hard rule is not formally broken. But "family = what the aligner could not separate" makes O2's
  within-family scope true by construction instead of by biology.
- **What IS usable (opt-in, not a replacement):** multimapper links as extra **candidate** edges that must still
  pass the exon-share test. The union `S|C` is precision-neutral on all three substrates. Against the same
  operator it gains +.033 (gorilla), +.008 (chr16) and 0 (dev), which is small. The larger recoverable class
  is **retrocopies** (mechanism 1). A targeted fix would score the exon-share test on the spliced transcript
  rather than the genomic span, which is the thing the reads get right. It is a DNA-side change and needs no
  reads.
- **One-line answer for the advisor:** *read Jaccard finds 2–7× more paralog pairs, but on segmental duplications
  it cannot tell a gene family from genes that were duplicated together (chr16 precision .96 → .72). The exon-sum
  fraction is the test that separates the two. Chaining the multimapper anchors to rebuild that test from reads
  recovers the precision but loses the recall, because reads never cover a gene's whole exon set.*

## §5 Caveats

- Guided nodes only (annotated genes). In de novo mode the loci come from reads, so the expression ceiling in §2
  is already paid.
- One library per species (human A119b, gorilla testis). Jaccard values will change with depth and tissue; the
  DNA rule's will not.
- `S_shipped` vs read arms mixes operators (MCL vs components). The `S.cc` row is the like-for-like baseline and
  every conclusion above holds against it (gorilla gain larger, chr16 precision loss the same).
- The gorilla aggregate is dominated by PF43 (280-gene KRAB-ZNF). Without it: `S_shipped` .293 / .988,
  `J0.01` .407 / .786, `Cany` .365 / .982. The direction is unchanged.
- Proposed register row (not yet added; the register has uncommitted edits in another session): *"Read-set
  Jaccard of multimapping molecules (+ multimapper chaining anchors) replacing the exon-sum edge rule | 2–7×
  more truth pairs linked; gorilla F .208→.455 | chr16 (LCR16) family precision .964→.725: co-duplicated
  passengers (EIF3CL~NPIPB9, BOLA2/SLX1A block) share reads; anchor-chain exon coverage ≥.60 F .178/.128 <
  shipped | exon-share fraction is the co-duplication filter; reads usable only as candidate edges."*

## Reproduce

```bash
O=/mnt/linuxdisk/tmp/advisor_jaccard   # outputs; reads/anchors cached as <tag>.reads.pkl / <tag>.anchors2.pkl
python3 bench/mechanism/jaccard_anchor_test.py hsa20 /mnt/linuxdisk/tmp/ppar/chr20.bam /mnt/linuxdisk/tmp/sedef/chr20.genes.gff \
  /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa chr20 /mnt/linuxdisk/tmp/sedef/truth/chr20.tsv $O
python3 bench/mechanism/jaccard_anchor_test.py ggo44 /mnt/linuxdisk/tmp/gw22/sec/ggo44.bam /mnt/linuxdisk/tmp/gw22/sec/ref/NC_073244.2.genes.gff \
  /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta NC_073244.2 /mnt/linuxdisk/tmp/gw22/sec/ref/NC_073244.2.tsv $O
samtools view -b -o $O/chr16.bam /mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam chr16 && samtools index $O/chr16.bam
python3 bench/mechanism/jaccard_anchor_test.py hsa16 $O/chr16.bam /mnt/linuxdisk/tmp/sedef/chr16.genes.gff \
  /mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa chr16 /mnt/linuxdisk/tmp/sedef/truth/chr16.tsv $O
```
Per tag: `summary.tsv` (every arm), `ceiling.txt`, `pairs.tsv` (per pair: shared molecules any/multimap,
Jaccard, anchor-chain coverage, DNA edge, truth).
