# Pre-registration — per-copy NPIP and TBC1D3 recovery: our assembler vs the StringTie, FLAIR and isoseq collapse baselines, on the same reads

**Written 2026-09-30 00:05 PDT (the study was opened 2026-09-29 22:40), before any transcript model of any tool and
before any read at a family copy is opened.** What HAS been opened before this freeze: the annotations (CHM13 RefSeq
full GFF chr16/chr17, the two existing gorilla truths), the baselines' recipes and logs (which library, which flags),
our runs' keys/params/driver logs, one truth-vs-truth smoke test of the scorer (no tool model), and one read-support
smoke test on GAPDH (chr12, not a family copy). Two species, **human A119b (CHM13 v2.0) and gorilla OR6737 (GGO)**, are
reported separately and never pooled. StringTie, FLAIR and isoseq collapse are trusted baselines; the question is
where each method (ours included) loses a copy, not who wins. Nothing in `src/`, `bench/soto/` or `tools/` is edited,
no default is changed, nothing is committed. Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools/`.

## 1. Question

For every annotated copy of the two core-duplicon families (NPIP, TBC1D3), on the same alignments, which of the four
methods emits (a) any model at the copy, (b) a model with the copy's annotated intron chain, (c) a model that belongs
to that copy alone, and — for copies with reads but no model — at which step each method loses it. Genome-wide the
answer is known in aggregate (r1058–r1066, §6z8: on reference transcripts with ≥ 2 exact-chain reads we are the most
sensitive arm; with 1 read isoseq leads; the tools add single-read isoforms). It has never been laid out **per copy of
a multi-copy family**, which is the thesis's object.

## 2. Inputs and arms (all confirmed to be the same libraries)

| arm | file | library / recipe (confirmed from the recipe beside the output) |
|---|---|---|
| ours, default | `rustle_figures/runs/human_A119b/human_A119b.gtf` (252,217 tx, 09-25); `runs/gorilla_OR6737/gorilla_OR6737.gtf` (74,045 tx, 09-25) | streaming `copy_assign --assemble-only --genome-wide` on `A119b.t2t.bam` / `GGO_mm.bam` (`assemble.key`: args `[]`, bins `as_table 12dc766b` + `copy_assign 753a3b4d`). Effective settings from the run log and `params.tsv`: junctions **strict**, polish full (fraction 0.02, mono shadow, mono quantile 0.82, ISM ratio 0.7, retained ratio 10), pass-1 floor 2, coordinate dedupe on, seeding = primaries + secondaries at ≥ 0.98 × genome-wide best AS ("GOOD" secondaries, §6z7). This is the shipped default as of commit 3007c3d4 (`--assemble-only` = strict + retained-ratio 10 since r1074–r1076); the working tree's flipped defaults are NOT used. `ggo_npip_loss` and `ggo_tbc1d3_holdout` showed the frozen `cc_bin_frozen/copy_assign` (57168f3f) reproduces this GTF byte for byte on the five gorilla family contigs; no such check exists for human (disclosed). |
| ours, primaries only | `human_A119b.primary.gtf` (238,406 tx); `gorilla_OR6737.primary.gtf` (72,690 tx) | the same with `--no-seed-secondaries` (`assemble_primary.key`) |
| StringTie 3.0.1 | `benchmark_collapse/stringtie_A119b/A119b.stringtie.gtf` (250,500 tx); `stringtie_GGO/GGO.stringtie.gtf` | `stringtie -L -p 8` on the FULL `A119b.t2t.bam` / `GGO_mm.bam`, no `-G` (`run_stringtie.sbatch`; log confirms bam=A119b.t2t.bam). Secondaries never seed an isoform (recipe note); 1-read floor for multi-exon isoforms. |
| FLAIR 3.0.1 | `flair_A119b/A119b.flair.isoforms.gtf` (1,096,281 isoforms); `flair_GGO/GGO.flair.isoforms.gtf` (153,402) | `flair_bam2bed.py` = FLAIR's `dofiltering()` on the same BAMs (primaries only, supplementary removed, `--quality 0`), then `flair collapse --trust_ends --generate_map`, no correct, no annotation (`run_flair.sbatch`). `--trust_ends` forces quality 0, so **no MAPQ filter**; the collapse support floor is FLAIR's default (`-s 3`). Per-isoform read counts in `*.isoform.counts.txt`. |
| isoseq collapse | `isoseq_upload/isoseq_A119b/A119b.collapsed.sorted.gff.gz`; `isoseq_GGO_OR6737/GGO_OR6737.collapsed.sorted.gff.gz` | `isoseq collapse` on the same mapped BAMs (`isoseq_collapse_*.sbatch`: MAPPED=A119b.t2t.bam / GGO_mm.bam), defaults `--min-aln-coverage 0.99 --min-aln-identity 0.95 --max-fuzzy-junction 5 --max-5p-diff 50 --max-3p-diff 100`; secondaries discarded on 0x100 (recipe note); no read floor. FL counts in `*.abundance.txt.gz`. |
| reads | `A119b.t2t.bam` (96.3 GB), `GGO_mm.bam` (11.7 GB) | minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`; primaries = `-F 2308` |

r867's trap (tool numbers from a different library) does not apply: every arm above ran on the same two BAMs.

## 3. Truths — annotation only (`code/build_truth.py` 5b2a6d82; `ann/truth.{hsa,ggo}.json`, `.gtf`, `copies.*.tsv`)

Construction mirrors `ggo_npip_sim/build_tx.py` (FL): territory = merged union of the record's raw annotated exons
(0-based half-open); per transcript the exons are canonicalised (every non-canonical annotated intron is merged;
canonical = GT-AG, GC-AG, AT-AC on the transcript strand — §6m7's rule: filter truth to canonical motifs first);
junction = (last exon base, first exon base) 1-based; primary chain = most junctions, then longest, then id; the truth
GTF carries one transcript per distinct canonical chain (terminal exons run to the territory ends, which gffcompare's
`=` ignores for multi-exon chains). Every chain of every copy is multi-exon (0 single-exon truth chains).

**Human NPIP — 26 copies** = the Dishuck-checked chr16 copies of `docs/lit_subclusters_npip_dishuck_check.tsv`
(NPIPB1P on chr18 left out), matched by GeneID (the 22 NPIP rows of `docs/lit_subclusters_npip_tbc1d3_truth.tsv` are a
subset: no PKD1P6-NPIPP1, LOC128966608, LOC124907834, LOC124907808, LOC124907807). Two readthrough-defined copies,
flagged `readthrough=1` — their only annotated isoforms span an SD partner, so a COMPLETE model there is FUSED by
construction: **PKD1P6-NPIPP1** (transcribed pseudogene = the readthrough; territory = its NPIPP1 half, exons at or
above chr16:15,126,650; isoforms = its own 2 transcripts, 28-junction primary) and **NPIPB14P** (a pseudogene record
with no exon features; territory = the exons of PDXDC2P-NPIPB14P inside its span; isoforms = that readthrough's
transcript, 25 junctions). Biotypes: 19 protein-coding, 2 pseudogene (NPIPB10P, NPIPB14P), 5 transcribed pseudogene
(PKD1P6-NPIPP1, LOC128966608, LOC124907834, and none else — LOC124907808/807 are protein-coding). 147 human truth
transcripts in all; non-canonical annotated introns merged: NPIPB13 11, NPIPA2 5, NPIPB4 5, NPIPB7 4, NPIPA5 4, others ≤ 3.

**Human TBC1D3 — 16 chr17 records** = every chr17 gene/pseudogene whose description matches "TBC1 domain family member
3" in `Reference/chm13v2.0_RefSeq_full.gff.gz` (never `HSA_genomic.gff`): the 9 protein-coding copies TBC1D3, 3B, 3D,
3E, 3F, 3G, 3H, 3I, 3K (13–14-junction primaries; = the 9 TBC1D3 rows of the lit truth), the 2 transcribed
pseudogenes TBC1D3P5 (12 junctions) and TBC1D3P2 (13), and **5 pseudogene records with no exon features** (TBC1D3P4,
P3, P7, P1, LOC100420311): their territory is the record span; they have no annotated isoform, so COMPLETE is **not
scorable** there (`n/s`), and they count only for LOCUS and read support.

**Gorilla NPIP — 25 copies** = the T_member truth `ggo_npip_sim/ann/truth.json`, reused verbatim (territory = native
exons inside the landing ∪ lifted exons; chains = the native record's transcripts clipped to the territory; 9 copies
have clipped chains — NPIPA1, NPIPB1P, NPIPB4, NPIPB8, LOC124907807, NPIPA5, NPIPA6, LOC124907834, NPIPB15 — so a
model reproducing the full native transcript there scores `k`/`j`, not `=`; reported). **Gorilla TBC1D3 — 14 records**
= `ggo_tbc1d3_holdout/ann/truth.json` (`PREREG_ggo_tbc1d3_holdout_2026-09-29.md` §2: 12 on NC_073228.2, the 3B-like
pseudogene + 3K-like pair on NC_073224.2), reused verbatim. 146 gorilla truth transcripts.

Contigs scored: human chr16 + chr17; gorilla NC_073241.2, NC_073242.2, NC_073244.2 (NPIP), NC_073228.2, NC_073224.2
(TBC1D3). Every tool's set is restricted to these contigs before gffcompare (`code/gc.sh` 5780c009, `gffcompare
v0.12.10 -r truth.gtf`, the house call of `figures/assembly.py` / `bench/score.py members`).

## 4. Definitions (per copy, per tool; `code/models.py` 724be65f, `code/reads.py` d0c63b99)

- **(a) LOCUS**: ≥ 1 model with ≥ 1 bp exonic overlap with the copy territory on the copy strand. Unstranded models
  (strand `.`, StringTie single-exon) count and are flagged. (The `npf_audit` presence definition.)
- **(b) COMPLETE**: ≥ 1 LOCUS model whose gffcompare class is `=` against a truth transcript **of this copy** (tmap
  `ref_gene_id` = the copy id). Cross-check `complete_exact`: the model's intron chain equals one of the copy's
  canonical chains (Python, no gffcompare); every disagreement is listed. Also **partial-or-better** = best class in
  {`=`, `c`, `k`} (r1062's "partial" set), and the best class per copy.
- **(c) UNIQUE / FUSED**: FUSED = ≥ 1 LOCUS model whose exons also overlap (same strand, ≥ 1 bp) another truth copy's
  territory, or the exons of any annotated RefSeq record (`families_gw/species/{human,gorilla}/genes.tsv`) that lie
  **outside** the copy territory span [lo, hi) — a partner. Records nested inside the span never make a partner.
  UNIQUE = LOCUS and not FUSED. The number of FUSED models per copy and the partners are reported.
- **models per copy** = LOCUS models (over-splitting), median and max over LOCUS copies, per tool.
- **(d) read support** (BAM, primaries `-F 2308`): own primary = orientation (FLAG 0x10 under `-uf`) equals the copy
  strand and ≥ 1 aligned base on the territory (the `ggo_npip_loss` definition). Per copy: own primaries; of which
  MAPQ ≥ 1; spliced; **exact-chain** = own primaries whose intron chain equals a canonical annotated chain of the copy
  (and after the shipped coordinate dedupe); `max_same_chain` = the largest number of own primaries sharing one intron
  chain (the grouping view of FLAIR's `-s 3` and isoseq's floor 1); whole-chain (junction set ⊇ primary chain, the
  trace's definition); isoseq-eligible (identity ≥ 0.95 and query coverage ≥ 0.99, approximating collapse's filter);
  other-orientation primaries; own secondaries (information).
- **(e) first losing step**: for OURS, the `ggo_npip_loss` port — gorilla NPIP and TBC1D3 rows are reused verbatim from
  `ggo_npip_loss/score/table.md` and `ggo_tbc1d3_holdout/diag/score/table.md` (OR rows; steps (0) no reads, (1)
  seeding, (2) dedupe, (4) pass-1 floor, gate (3) strict junctions / (8) single-exon `+` placeholder, (6) polish mono
  floor, (7) merged); human copies absent from our GTF are traced with the same port (`code/trace_port.py` 25bcd428, a
  path-only copy of the frozen `trace.py`: seeding → dedupe → pass-1 → gate; a copy with gate survivors but no shipped
  model is "polish" and, if any occurs, the no-polish regional arm is run to name the polish step). For the BASELINES,
  what their inputs/outputs expose at a copy they miss: own primaries (MAPQ ≥ 1, spliced), `max_same_chain` (vs FLAIR's
  floor 3), isoseq-eligible reads (vs collapse's filter), and for present copies the model's own support (StringTie
  `cov`, FLAIR read count, isoseq `count_fl`, ours `reads`).

**Denominators** (fixed from the BAM before any model is opened): ALL = every copy; **R2 = copies with ≥ 2 own
primaries** (the headline); R0 = copies with 0 own primaries; for COMPLETE: scorable copies (with ≥ 1 annotated
isoform), R2 ∩ scorable, **E2 = copies with ≥ 2 exact-chain own primaries** (r1065's universe), E1 = exactly 1. No
denominator is conditioned on any prediction (metric trap 8). Sensitivity = recovered / denominator, reported as
counts and fractions; no test statistics at n ≤ 26.

## 5. Predictions (numeric, per species × family cell; ours = the default arm)

Prior read-level knowledge, disclosed: gorilla OR6737 NPIP (ours) 11/25 LOCUS, 0 COMPLETE, per-copy steps known
(`ggo_npip_loss`); gorilla OR6737 TBC1D3 (ours) 7/14 LOCUS, 6 COMPLETE, two present copies (LOC129533806,
LOC129533813) have 0 own primaries, three absent ones have 0 primaries, four absent ones have 1–6 primaries each in a
distinct chain (`ggo_tbc1d3_holdout`); human chr16 NPIP at a seeded BASE run 26/26 LOCUS with 11 copies in fused loci
(`npf_audit`); genome-wide r1065. **Nothing is known about any baseline's models at these copies, nor about human
TBC1D3 for any method.**

- **P1 (floors → LOCUS).** In every cell, isoseq collapse's LOCUS count on R2 ≥ ours; StringTie's ≥ ours on gorilla
  NPIP (our absences there are pass-1-floor/gate losses of copies with ≥ 2 primaries); FLAIR's ≤ isoseq's everywhere.
- **P2 (seeding).** Every baseline has LOCUS at 0 copies of R0 in every cell (they are primary-only); ours has LOCUS at
  ≥ 1 R0 copy (gorilla TBC1D3's two). Ours-primaries-only ⊆ ours in LOCUS count in every cell.
- **P3 (r1062–r1066 per copy).** On E2, our COMPLETE count ≥ every baseline's in every cell where |E2| ≥ 1 (ties
  allowed). Expected testable cells: human NPIP, human TBC1D3, gorilla TBC1D3; gorilla NPIP's E2 is expected empty.
- **P4 (single exact read).** On E1, ours ≤ isoseq in every cell; isoseq ≥ 1 in any cell with |E1| ≥ 2.
- **P5 (gorilla NPIP completeness is lost before the assembler).** COMPLETE = 0 for **every** method at gorilla NPIP
  (OR6737 has 0 whole-chain reads at all 25 copies).
- **P6 (fusion is in the reads).** Human NPIP: every method, ours included, has FUSED at ≥ 3 of the 4 RefSeq-annotated
  readthrough copies NPIPA1 (PKD1P3-NPIPA1), NPIPA6 (PKD1P1-NPIPA5L), NPIPA9 (PKD1P5 readthrough), NPIPB14P
  (PDXDC2P-NPIPB14P); ours 8–12 FUSED copies in all. Gorilla: ≤ 1 FUSED copy per method per family.
- **P7 (over-splitting).** Median models per LOCUS copy: isoseq ≥ ours in both human cells; FLAIR ≤ isoseq.
- **P8 (human NPIP presence).** Ours 26/26 LOCUS; every baseline ≥ 24/26.
- **P9 (what the baselines' misses look like).** At copies with ≥ 1 own primary that FLAIR misses, `max_same_chain`
  < 3 in ≥ 80% of the cases; at copies isoseq misses, the miss is not explained by the eligibility filter alone in ≥ 50%
  of cases (isoseq-eligible reads ≥ 1) — i.e. isoseq's per-copy losses will be placement, not filtering.

## 6. Falsifiers (each is filed as a register row if it fires)

- **F1**: a cell where a baseline's COMPLETE on E2 exceeds ours by ≥ 2 copies, or by ≥ 1 copy in ≥ 2 cells ⇒ the
  aggregate r1065 claim does not hold per copy on the thesis's families.
- **F2**: any `=` model at any gorilla NPIP copy by any method ⇒ `ggo_npip_loss`'s "lost before the assembler" is
  refuted or its whole-chain instrument is wrong; the read(s) behind the model are named.
- **F3**: a stranded baseline model with LOCUS at an R0 copy ⇒ the primary-only claim for that tool, or my own-read
  definition, is wrong; which one is stated.
- **F4**: ours LOCUS on R2 below every baseline by ≥ 3 copies in any cell ⇒ the floor/gate precision trade costs more
  at the thesis's families than the genome-wide numbers show.
- **F5**: gffcompare `=` and `complete_exact` disagree on > 2 models in any cell ⇒ scorer defect, fixed as an
  instrument and reported before any verdict.

## 7. Hostile self-review (written before scoring)

1. **The truth is "in width" only.** Copies absent from the annotation are invisible; a method's models at unannotated
   copies are neither credited nor penalised. Nothing here is O3.
2. **`=` is exact on the intron chain.** A model whose junctions sit 1–5 bp off (isoseq's fuzzy-junction merging keeps
   a read's own junctions, so this is rare) scores `j`; the best-class column and partial-or-better show it.
3. **LOCUS at ≥ 1 bp overlap is lenient**: a 3′ single-exon fragment counts for every method (the ape NPIP loci are
   mostly such fragments, `ggo_npip_loss`). LOCUS is a presence bound, not a quality claim; COMPLETE is the quality claim.
4. **Two human copies are readthrough-defined** and can only be COMPLETE and FUSED at once; they are flagged and the
   FUSED counts are reported with and without them.
5. **Clipped gorilla NPIP chains** (9 copies) can turn a correct full-length model into `k`; P5 predicts 0 regardless,
   and `k` is reported.
6. **Unstranded StringTie models** are counted for LOCUS on either strand — a small lenience for StringTie only; the
   count of unstranded models is reported.
7. **The exact-chain read count is computed from the same BAM every tool saw**, so E2/E1 are method-independent; but
   FLAIR and isoseq also see the reads' 5′/3′ ends, which the count ignores — a copy in E2 can still be split by ends.
8. **Different tissues and depths across species** (A119b 21.4 M molecules; OR6737 4.4 M): species are never pooled
   and no cross-species sensitivity is compared.
9. **The 0.98-AS seeding gives ours reads the baselines never see** (P2). Where ours is present at an R0 copy the
   presence rests on reads whose best placement is elsewhere; that is reported as such, not as a recovery.
10. **The five no-exon TBC1D3P records** are scored for LOCUS only; a method's model there is presence at a pseudogene
    span, nothing more.
11. **The first-losing-step port for human** reuses the gorilla instrument unchanged (path-only copy); it was validated
    against the binary on gorilla (I2/I3 of the sibling studies), not on human. If a human copy with reads is absent
    and the port's gate survivors disagree with the shipped GTF, the regional frozen-binary arms are run before the step
    is named.
12. **n is small** (14–26 copies per cell). Counts are reported; no fraction is called a difference without the counts.

## 8. Frozen code, stop rules

`code/` sha1 (first 8): `build_truth.py` 5b2a6d82 · `reads.py` d0c63b99 · `models.py` 724be65f · `gc.sh` 5780c009 ·
`report.py` 666b400b · `trace_port.py` 25bcd428 (diff vs `ggo_npip_loss/trace.py`: the five path lines only). Machine:
light lock for every gffcompare / python step; the heavy lock only if a regional frozen-binary arm is needed (§4e).
**After this freeze nothing in §3–§6 changes.** Instrument defects are fixed, recorded in the Outcome with the sha1
before/after, and never used to move a bar. Register rows are drafted with suffix G and not appended. Outcome below.

## Outcome

*(Appended 2026-09-30 00:19 PDT. The text above this heading is byte-identical to the frozen version, sha1 40b581b8,
recorded in `FROZEN.sha1` at 00:06 before any tool model or family read was opened; re-verified before this append.
Everything was run as registered; the three post-freeze instruments are listed in O6.)*

### O1. Answer first

- **P3 holds per copy — where ≥ 2 primary reads carry an annotated chain, ours is COMPLETE at least as often as every
  baseline in all three testable cells** (human NPIP 11/13 vs StringTie 10, isoseq 9, FLAIR 6; human TBC1D3 9/9 vs
  9/9/8; gorilla TBC1D3 3/3 vs 3/3/2). Gorilla NPIP has 0 such copies and **0 COMPLETE for every method (P5)**.
- **F4 fired at gorilla NPIP**: LOCUS on the 23 copies with ≥ 2 primaries is ours 11 vs StringTie 16, FLAIR 18,
  isoseq 23. Qualifier (post-hoc, O6): a model carrying ≥ 1 annotated junction exists at StringTie 3/16, FLAIR 3/19,
  isoseq 10/24 of their LOCUS copies (ours 0/11); isoseq's best model is single-exon at 19 of its 24. The extra
  presence is fragments; no method reconstructs a gorilla NPIP isoform.
- **P4's second clause failed**: the single-exact-read chains are recovered by **StringTie** (its 1-read floor; `=`
  at cov 1–2 at human NPIPB12, NPIPB13 and gorilla LOC129533792), not by isoseq (0). FLAIR also matches gorilla
  LOC129533792 at "support 3" by folding two 5′-truncated sub-chain reads (10 and 5 of 13 junctions) into the one
  full-chain read — the 5′ fold-in that r1071 refuted for us.
- **Our two human NPIP E2 misses are polish steps, named with the no-polish regional arm** (I1 byte-identical to the
  09-25 genome-wide GTF on chr16+chr17; polish port ALL IDENTICAL): NPIPB2's two truth-chain transcripts (4 and 2
  reads) fall to the isoform-fraction rule because the locus is fused with GSPT1 (948 reads; r1063's family);
  NPIPB6's 2-read truth chain is an ISM sub-chain of a 67-read `j` model (r1064's family).
- **Seeding is what separates us from the baselines at low-read copies**: 6 of our COMPLETE copies have < 2
  exact-chain primaries (gorilla LOC129533808/806/813, LOC129533797; human TBC1D3B, TBC1D3I) and rest on
  GOOD secondaries; every baseline (primary-only) has no locus at the two 0-primary gorilla copies (P2 ✓) and only
  `j`/`n`/`m` fragments at TBC1D3B/I. The primaries-only arm loses exactly these (human TBC1D3 11→9, gorilla 7→3).
- **Fusion is in the reads (P6)**: all five arms are FUSED at the four RefSeq readthrough copies; human NPIP FUSED
  copies ours 13 (11 without the two readthrough-defined copies; the bar said 8–12), StringTie 14, FLAIR 17,
  isoseq 20. The gorilla clause failed: isoseq is FUSED at 11 gorilla NPIP copies (fragments joined to neighbour
  LOCs), StringTie 4, FLAIR 3, ours 1.
- **Over-splitting (P7 ✓)**: median models per present copy, human NPIP: ours 16 (max 76), primaries-only 7,
  StringTie 13, FLAIR 20, isoseq 74 (max 321); human TBC1D3: 12 / 6 / 12 / 14 / 44.

### O2. Per-family per-tool sensitivities (counts; denominators fixed from the BAM)

Human A119b — **NPIP, 26 copies** (all 26 have ≥ 2 own primaries; 13 have ≥ 2 exact-chain reads, 2 exactly 1, 11 none):

| arm | LOCUS | COMPLETE /26 | COMPLETE on E2 /13 | on E1 /2 | partial-or-better (=,c,k) | FUSED (excl. RT-defined) | UNIQUE | models per LOCUS copy med/max |
|---|---|---|---|---|---|---|---|---|
| ours (default) | 26 | **11** | **11** | 0 | 15 | 13 (11) | 13 | 16 / 76 |
| ours, primaries only | 26 | 10 | 10 | 0 | 12 | 12 (10) | 14 | 7 / 42 |
| StringTie 3.0.1 -L | 26 | 12 | 10 | **2** | 15 | 14 (12) | 12 | 13 / 37 |
| FLAIR 3.0.1 | 26 | 6 | 6 | 0 | 12 | 17 (15) | 9 | 20 / 92 |
| isoseq collapse | 26 | 9 | 9 | 0 | 21 | 20 (18) | 6 | 74 / 321 |

Human A119b — **TBC1D3, 16 records** (all ≥ 2 primaries; 11 scorable for COMPLETE, 9 with ≥ 2 exact-chain reads, 0 with 1):

| arm | LOCUS /16 | COMPLETE /11 | on E2 /9 | partial-or-better | FUSED | UNIQUE | models med/max |
|---|---|---|---|---|---|---|---|
| ours (default) | 13 | **11** | 9 | 11 | 5 | 8 | 12 / 21 |
| ours, primaries only | 13 | 9 | 9 | 9 | 6 | 7 | 6 / 20 |
| StringTie | 14 | 9 | 9 | 9 | 5 | 9 | 12 / 41 |
| FLAIR | 13 | 9 | 9 | 9 | 9 | 4 | 14 / 39 |
| isoseq | 16 | 8 | 8 | 9 | 8 | 8 | 44 / 166 |

Gorilla OR6737 — **NPIP, 25 copies** (24 with ≥ 1 primary, 23 with ≥ 2, NPIPA7 has 0; **0 copies with an exact-chain read**):

| arm | LOCUS /25 | LOCUS on R2 /23 | COMPLETE | partial-or-better /23 | any model with an annotated junction (post-hoc) | best model multi-exon | FUSED | models med/max |
|---|---|---|---|---|---|---|---|---|
| ours (default) | 11 | 11 | 0 | 1 | 0 / 11 | 10 / 11 | 1 | 2 / 3 |
| ours, primaries only | 9 | 9 | 0 | 1 | 0 / 9 | 8 / 9 | 1 | 1 / 1 |
| StringTie | 16 | 16 | 0 | 6 | 3 / 16 | 10 / 16 | 4 | 3 / 10 |
| FLAIR | 19 | 18 | 0 | 7 | 3 / 19 | 11 / 19 | 3 | 2 / 7 |
| isoseq | 24 | 23 | 0 | 19 | 10 / 24 | 5 / 24 | 11 | 4 / 21 |

Gorilla OR6737 — **TBC1D3, 14 records** (9 with ≥ 1 primary, 7 with ≥ 2, 5 with 0; 3 with ≥ 2 exact-chain reads, 1 with exactly 1):

| arm | LOCUS /14 | on R2 /7 | at 0-primary copies /5 | COMPLETE /14 | on E2 /3 | on E1 /1 | FUSED | models med/max |
|---|---|---|---|---|---|---|---|---|
| ours (default) | 7 | 4 | **2** | **7** | 3 | 0 | 1 | 1 / 7 |
| ours, primaries only | 4 | 4 | 0 | 3 | 3 | 0 | 1 | 1 / 6 |
| StringTie | 7 | 7 | 0 | 4 | 3 | 1 | 1 | 2 / 6 |
| FLAIR | 5 | 5 | 0 | 4 | 3 | 1 | 1 | 1 / 5 |
| isoseq | 9 | 7 | 0 | 2 | 2 | 0 | 2 | 3 / 27 |

gffcompare `=` and the Python exact-chain check agree on every model in every cell (0 disagreements; F5 not fired).
No arm emits an unstranded model at any copy. The full per-copy tables (read support, per-tool L/C/F/n, best class,
support, our first losing step) are in the report `scratchpad/figs/copy_recovery_tools.md` and `score/report.md`.

### O3. Predictions and falsifiers

| item | verdict | numbers |
|---|---|---|
| P1 floors → LOCUS | **✓** all clauses | isoseq ≥ ours on R2 in 4/4 cells (26=26, 16>13, 23>11, 7>4); StringTie 16 > ours 11 at gorilla NPIP; FLAIR ≤ isoseq in 4/4 |
| P2 seeding | **✓** | baselines 0 LOCUS at R0 copies in every cell (gorilla NPIP 0/1, TBC1D3 0/5); ours 2/5 at gorilla TBC1D3; primaries-only ⊆ default in LOCUS in 4/4 |
| P3 E2 completeness | **✓** | ours ≥ every baseline in the 3 cells with \|E2\| ≥ 1; gorilla NPIP E2 = ∅ as predicted |
| P4 E1 | **half** | ours ≤ isoseq (0 ≤ 0) ✓; "isoseq ≥ 1 when \|E1\| ≥ 2" **✗** (0/2 at human NPIP — StringTie 2/2 instead) |
| P5 gorilla NPIP | **✓** | 0 COMPLETE for all five arms |
| P6 fusion | **partly** | FUSED at ≥ 3 of the 4 readthrough copies: 4/4 for every arm ✓; ours 8–12 ✗ (13; 11 excl. RT-defined); gorilla ≤ 1 per method ✗ (isoseq 11 / 2, StringTie 4 / 1, FLAIR 3 / 1; ours 1 / 1) |
| P7 over-splitting | **✓** | isoseq median 74 vs ours 16 (NPIP), 44 vs 12 (TBC1D3); FLAIR ≤ isoseq |
| P8 human NPIP presence | **✓** | 26/26 for every arm |
| P9 baseline misses | **✗ / vacuous** | FLAIR misses at R1 with max-same-chain < 3: 8/12 = 67% (< 80%; NPIPA1 22 primaries/15 same-chain and NPIPB11 28/8 are FLAIR misses its outputs do not explain); isoseq misses 0 copies with ≥ 1 primary in every cell |
| F1 | not fired | no cell where a baseline beats ours on E2 |
| F2 | not fired | 0 `=` at gorilla NPIP |
| F3 | not fired | no baseline model at any 0-primary copy |
| **F4** | **FIRED (gorilla NPIP)** | ours 11 vs 16 / 18 / 23 on R2 (−5 to −12); the other three cells do not fire (human TBC1D3 13 vs 14/13/16; gorilla TBC1D3 4 vs 7/5/7 — FLAIR only +1) |
| F5 | not fired | 0 gffcompare/exact disagreements |

### O4. Where each method loses a copy

- **Ours, gorilla NPIP (14 absences, steps reused from `ggo_npip_loss`)**: pass-1 floor 3 (NPIPB10P, NPIPA7, NPIPA5),
  gate 8 (NPIPB1P, NPIPA8, NPIPB8, NPIPA9, LOC124907808, NPIPB2, NPIPB6, NPIPB15), mono floor 3 (NPIPB14P, NPIPB12,
  NPIPB5). At these 14, StringTie has a locus at 6, FLAIR at 10, isoseq at 13 (support 1–25 reads); a model with an
  annotated junction exists there only at NPIPB14P (StringTie, FLAIR, isoseq), NPIPA5 (FLAIR, isoseq), NPIPB5 and
  NPIPB15 (isoseq). NPIPB8 (66 primaries, lost at strict junctions) is `c` for all three baselines with 9–25 reads
  and 0 annotated junctions.
- **Ours, gorilla TBC1D3 (7 absences, reused)**: pass-1 floor 4 (1–6 primaries, every read a distinct chain: the
  baselines are present there with 1-read (StringTie/isoseq) or fold-in (FLAIR) models), seeding 3 (0 primaries;
  every baseline absent too).
- **Ours, human (3 absences, all TBC1D3 pseudogene spans; trace port + no-polish arm)**: TBC1D3P4 (2 primaries): the
  only gate survivor on the span is a 593-read single-exon cluster dropped as a mono **shadow**; TBC1D3P3 (4): gate
  (strict junctions + single-exon `+` placeholder on a `−` span); TBC1D3P7 (2): 2-read single-exon cluster below the
  mono **floor** (chr17 floor 11). StringTie misses P3 and P7, FLAIR misses all three, isoseq has 2–3 one-read
  single-exon models at each.
- **Ours, incomplete at E2 (2, human NPIP)**: NPIPB2 — truth-chain transcripts of 4 and 2 reads dropped by
  **fraction** (0.02 of a 948-read locus fused with GSPT1; the copy's own reads are 23 of it); NPIPB6 — the 2-read
  truth chain is a sub-chain of a 67-read model, dropped by **ISM** (ratio 0.7). Both steps are closed rows
  (r1063, r1064); the NPIPB2 loss is the fusion cost `npf_audit` described (GSPT1 ghost bridge).
- **Baselines at copies with reads but no model**: FLAIR — 12 misses across cells, 8 with max-same-chain < 3 (its
  floor); the 4 others (NPIPA1, NPIPB11, NPIPA2, NPIPB6 in gorilla) are not explained by anything its outputs
  expose. StringTie — 10 misses, 9 of them copies with ≤ 9 primaries and max-same-chain ≤ 5 (NPIPA1, 22 primaries /
  15 same-chain, is the exception). isoseq — misses no copy that has a primary read (P9 second clause vacuous); its
  loss is completeness, not presence.

### O5. Hostile self-review, revisited

1. LOCUS ≥ 1 bp is lenient, as written: at gorilla NPIP it credits isoseq with 24 loci of which 19 are single-exon
   fragments; the annotated-junction qualifier (O2) is the honest presence measure and was added post hoc (O6).
2. `=` is exact: the gorilla NPIP clipped chains (9 copies) did not matter — no method gets `=`, `k` or `j` with an
   annotated junction there except the four copies named in O4.
3. Two human copies COMPLETE only through seeded secondaries (TBC1D3B, TBC1D3I; 0 exact-chain primaries) count as
   COMPLETE by the pre-registered definition and are reported as resting on reads placed better elsewhere, per §7.9.
4. The E1/E2 split is from primaries; FLAIR's LOC129533792 `=` shows a baseline can reach `=` from 1 exact read + 2
   sub-chain reads. The split is still method-independent; it just is not FLAIR's own support count.
5. n ≤ 26 per cell; every claim above is a count.

### O6. Deviations and post-freeze instruments (none changes a definition, denominator or bar)

1. `code/arm_hsa.sh` (regional frozen-binary `ship`/`nopol` on chr16+chr17, human paths of `ggo_npip_loss/arm.sh`;
   heavy lock, 52 s / 0.78 GB and 23 s / 0.75 GB). **I1: the regional `ship` GTF equals the 09-25 genome-wide GTF on
   chr16+chr17 line for line, cov/TPM stripped (201,183 lines).**
2. `code/polish_port.py` (path-only copy of the frozen port, one line): **ALL IDENTICAL** on chr16 (floor 10; kept
   9,473 = shipped) and chr17 (floor 11; 10,668 = shipped).
3. `code/posthoc.py`: the annotated-junction / multi-exon qualifiers, the P6/P9 counts, the FLAIR LOC129533792 read
   lookup, and the isoseq-only gorilla NPIP loci table (`score/posthoc.json`).
4. The prereg named §4e's human step from the trace port; for the two "polish" copies and the two E2 misses the
   no-polish arm was run as §4e/§7.11 required.

### O7. Register rows (filed 2026-09-30 in `docs/NEGATIVE_RESULTS_REGISTER.md` as 1189-1193; drafted as G1-G5)

| id | date | § | claim tested | outcome |
|---|---|---|---|---|
| 1189 | 2026-09-30 | copy recovery vs baselines | Per copy of NPIP/TBC1D3 on the same reads, our shipped floor 2 / strict junctions / mono floor cost no presence relative to StringTie, FLAIR and isoseq collapse (`PREREG_copy_recovery_tools_2026-09-29.md`, sha1 40b581b8, F4) | ⛔ **F4 fired at gorilla NPIP: LOCUS on the 23 copies with ≥ 2 primaries is ours 11 vs StringTie 16, FLAIR 18, isoseq 23.** But the extra loci are fragments: a model with ≥ 1 annotated junction exists at StringTie 3/16, FLAIR 3/19, isoseq 10/24 of their present copies (ours 0/11), isoseq's best model is single-exon at 19/24, and every method is 0 COMPLETE at all 25 copies. Human NPIP/TBC1D3 and gorilla TBC1D3 do not fire (13 vs 14/13/16; 4 vs 7/5/7 on R2). Nothing here reopens r1067–r1073 / r1074–r1076: presence there is bought with unspliced or unannotated fragments. |
| 1190 | 2026-09-30 | copy recovery vs baselines | isoseq collapse (no floor) recovers the annotated chain at copies whose ONLY exact-chain read is single (P4) | ⛔ **isoseq 0/2 (human NPIP E1) and 0/1 (gorilla TBC1D3 E1); StringTie's 1-read floor does: `=` at NPIPB12 (cov 1.0), NPIPB13 (2.0) and gorilla LOC129533792 (1.0); FLAIR reaches the gorilla one at "support 3" by folding two 5′-truncated sub-chain reads (10 and 5 of 13 annotated junctions) into the one full read (the 5′ fold-in r1071 refuted for us).** Ours 0 by construction (pass-1 floor 2). On E2 (≥ 2 exact-chain primaries) ours ≥ every baseline in all 3 testable cells (11/13, 9/9, 3/3 vs best baseline 10/13, 9/9, 3/3) — r1065's aggregate holds per copy at the thesis's families. |
| 1191 | 2026-09-30 | copy recovery vs baselines | Our two human NPIP copies with ≥ 2 exact-chain reads but no `=` model (NPIPB2, NPIPB6) are assembly losses | ⛔ **Both are polish steps, named with the no-polish regional arm (I1 byte-identical; polish port ALL IDENTICAL):** NPIPB2's truth-chain transcripts (4 and 2 reads) fall to the isoform-fraction 0.02 rule because its locus is fused with GSPT1 (948 reads, the copy's own 23) — the `npf_audit` ghost-bridge cost made concrete; NPIPB6's 2-read truth chain is an ISM sub-chain (ratio 0.7) of a 67-read `j` model. StringTie/isoseq are `=` at NPIPB6 and FLAIR/isoseq at NPIPB2 with 1–3 reads. Closed rows r1063 (fraction) and r1064 (ISM); not reopened. |
| 1192 | 2026-09-30 | copy recovery vs baselines | The baselines, being primary-only, cannot be present at copies whose reads all place better elsewhere; our GOOD-secondary seeding can (P2) | ⭐ **Confirmed: 0 baseline loci at any 0-primary copy (gorilla NPIP 0/1, TBC1D3 0/5); ours present AND complete at gorilla LOC129533806/813 (0 primaries) and complete at LOC129533808 (1 primary), human TBC1D3B/I (0 exact-chain primaries); the primaries-only arm loses exactly these (COMPLETE 7→3 gorilla TBC1D3, 11→9 human TBC1D3).** These rest on reads whose best placement is another copy (§7.9) and are reported as such, not as recoveries. |
| 1193 | 2026-09-30 | copy recovery vs baselines | Fusion at human NPIP is a defect of our locus formation rather than of the reads (P6) | ⭐/⛔ **All five arms are FUSED at the four RefSeq readthrough copies (NPIPA1, NPIPA6, NPIPA9, NPIPB14P); FUSED human NPIP copies ours 13 (11 excl. the two readthrough-defined copies) ≤ StringTie 14 < FLAIR 17 < isoseq 20; UNIQUE ours 13, StringTie 12, FLAIR 9, isoseq 6.** The gorilla clause failed the other way: isoseq is FUSED at 11 gorilla NPIP copies (one-read fragments joined to neighbour LOCs), StringTie 4, FLAIR 3, ours 1. Over-splitting: median models per present human NPIP copy isoseq 74 (max 321), FLAIR 20, ours 16, StringTie 13, ours-primaries 7. |

### O8. Status after the Outcome (2026-09-30)

- **Not independently recomputed.** Every count above comes from the session that ran the study. The `=` and
  exact-chain checks agree on every model (0 disagreements) but that is an internal check.
- **Recheck of the polish attribution (own code, GTFs and truth only):** in `runs/hsa.nopol.gtf` there are exact-chain
  models at NPIPB2 (2: `DN_chr16_11963325_10`, `DN_chr16_11963320_9`) and NPIPB6 (1: `DN_chr16_28623531_7`); in
  `runs/hsa.ours.gtf` (the shipped default) there are 0 at both. This confirms that both E2 misses are lost by a polish
  step.
- **Rows filed as 1189-1193** with two precision edits: G2 gained the all-26-copies tally (StringTie 12, ours 11,
  isoseq 9, FLAIR 6), and G4's "§7.9" became "prereg §7.9".
- **Consolidated write-up:** `docs/archive/2026-09/COPY_RECOVERY_TOOLS_2026-09-29.md`.
