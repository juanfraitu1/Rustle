# Locus formation (node width and composition) vs annotated gene records — agent 1 of 2

Declarations written before any number: `/mnt/linuxdisk/home/juanfraitu/locus_width/DECLARATIONS.txt`, 2026-09-18T02:15:26Z, md5 `ec5a4a2812186e266e6395b212614e9b`. Nothing in `src/` was modified, nothing committed. Outputs: `/mnt/linuxdisk/home/juanfraitu/locus_width/out/SUMMARY_locus.tsv` (every variant × substrate row), `out/locus_{a,b,c}.json`, `out/family_b.json`, scripts in `out/../scripts/`, new PAFs in `map/`.

> **Verifier corrections, applied by the orchestrator** (verifier ok = true, 11 corrections; the mechanism and the post-hoc read-isoform result survive):
> 1. **W4 does not reproduce on substrate (a)**: implementing the declared containment rule literally gives a different fold set than the builder's 2 folds.
> 2. **"W4 is never harmful" is too strong** — on the held-back gorilla substrate it costs 6 expressed records their node, and its locus-level gain is largely denominator shrinkage.
> 3. **The junction-presence column is not reproducible**: no declared convention, and six readings all disagree with the builder except on one cell. Do not quote `junc_present` absolutely; only its direction.
> 4. **"2664/2664 consolidate one piece" contradicts the same report's "4 wider"**: measured, 2,660 of 2,664 have exons == rep_exons on gorilla; 705 of 710 on the ideal substrate.
> 5. **NPIP denominator is 27, not 26**, on substrate (a) — all 27 copies are inside it.
> 6. **`nodes_no_rec` counts unmatched nodes, not nodes overlapping no record.** The informative quantity: nodes overlapping NO expressed annotated record are 31/710 (a), **282/560 (b)**, 464/2,664 (c).
> 7. **The W7 gain on real reads must be quoted with its merged-node count**: nodes overlapping ≥ 2 expressed records go 42 → 112 (k=2) / 102 (k=3); an independent reconstruction reproduces the effect (42 → 85).
> 8. **W6-W10 are post hoc and undeclared** — the declarations fix only W0-W5 and the annotated ceiling, so those arms are not reproducible from the declarations and must be pre-registered before being quoted.
> 9. Minor: the annotated ceiling arm has 15 tx keys without minimap2 hits (miRNAs of 57-79 bp, benign); 7 of 166 records differ in transcript-set extraction (metric 2 only).
> 10. **Trap watch**: FAMILY P NPIP-only is 1.000 in every row including the ones whose family structure collapses — vacuous here, as the metric-traps register says. Only P strict, F strict and FAMILY R discriminate.
> 11. Substrate (c) gorilla carries no weight in the decision rule, and is now burned for the post-hoc W9 arm because it was run there.

## 0. Parity, so the baseline is the shipped thing

- Substrate (b): my W0 node set is rebuilt from the run's own reps dump through the verified `shared_definition.rs` mirror and asserted **byte-for-byte equal, in order, to the frozen A4 node set** (534 reps → 522 gene-level loci + 38 read-locus nodes = 560 nodes). Its family scoring reproduces the frozen A4 row exactly: 1349 exon + 1885 body edges, 1919 pairs, 68 families / 364 loci, FAMILY R 1.000, P strict 0.290, F strict 0.450, 5/27 full-length.
- Substrate (a): W0 is the run's own node list (`nodes/nodes.IDEAL.tsv`, 710 nodes on chr16/17/18), with each node's consolidated pieces recovered the way `s12_fixmerged.py` recovered them for the merged nodes (every `[rep-audit]` rep inside the node span on the node strand, its chain resolved through `transcripts.tsv`). Node span == recovered exon-union span for 710/710. 47 nodes fall back to the table's exons because one of their reps could not be resolved (31 reps unresolved out of 697) — that is the one place where (a) is a reconstruction rather than a capture.
- Substrate (c): W0 is the shipped `nodes/GGO.dn.nodes.tsv` as written by the run (2664 nodes).

## 1. The mechanism, measured

The task statement says a node's `exons` is "the merged union of every piece consolidated into the locus". On all three substrates it almost never is:

| substrate | nodes | nodes with exactly 1 consolidated piece | nodes whose exon union is wider than the rep chain |
|---|---|---|---|
| (a) ideal | 710 | **705** | 3 |
| (b) testis | 560 | ~522 of 534 reps are 1:1 with a locus; W1 adds only **12** extra tx query keys over 560 | — |
| (c) gorilla | 2664 | **2664** | 4 |

The catalog emits one representative per locus *before* `shared_definition::consolidate` runs, so consolidation has nothing to merge. `exons` is therefore one isoform, the tx query is that same isoform, and the gene body is that isoform's span. That is why NPIP (4.0 transcripts per copy) loses full-length nodes while TBC1D3 (2.74) mostly does not. **W1 as declared — "the chains of every piece" — is a no-op in practice: 12 extra queries on (b), 0 new edges, every family metric bit-identical to W0.**

## 2. Declared variants: the full curve

Metric 1a is the declared primary (exon-bp Jaccard, one-to-one matching, evaluation only). W2 and W3 change the body only, so 1a is invariant under them by construction — this was declared in advance, which is why metric 1b (span form) was pre-registered alongside.

### (a) human ideal expression, chr16/17/18 — 710 nodes, 701 expressed records

| variant | nodes | matched | R | P | **F** | ΔF | pair recall | pair prec | tx contained | junctions |
|---|---|---|---|---|---|---|---|---|---|---|
| W0 | 710 | 673 | 0.960 | 0.948 | **0.954** | — | 0.858 | 0.991 | 0.603 | 0.698 |
| W1 | 710 | 673 | 0.960 | 0.948 | 0.954 | +0.000 | 0.858 | 0.991 | 0.603 | 0.698 |
| W2 | 710 | 673 | 0.960 | 0.948 | 0.954 | +0.000 | 0.858 | 0.991 | 0.603 | 0.698 |
| W3 E=0/250/500/1000/2000 | 710 | 673 | 0.960 | 0.948 | 0.954 | +0.000 | 0.858 | 0.991 | 0.603 | 0.698 |
| W4 | 708 | 673 | 0.960 | 0.951 | 0.955 | +0.001 | 0.858 | 0.991 | 0.603 | 0.698 |
| W5 | 708 | 673 | 0.960 | 0.951 | 0.955 | +0.001 | 0.858 | 0.991 | 0.603 | 0.698 |
| **ANN ceiling** | 701 | 701 | 1.000 | 1.000 | **1.000** | +0.046 | 1.000 | 1.000 | 1.000 | 0.772 |

Span form (1b): W0 0.955 → W3_E2000 0.958, W2 0.954. The extensions move it by ≤ 0.003.

### (b) human real testis, 18 windows — 560 nodes, 166 expressed records

| variant | nodes | matched | R | P | **F** | ΔF | pair recall | pair prec | spanF (1b) | tx contained | junctions |
|---|---|---|---|---|---|---|---|---|---|---|---|
| W0 | 560 | 139 | 0.837 | 0.248 | **0.383** | — | 0.569 | 0.566 | 0.402 | 0.132 | 0.338 |
| W1 | 560 | 139 | 0.837 | 0.248 | 0.383 | +0.000 | 0.569 | 0.566 | 0.402 | 0.132 | 0.338 |
| W2 | 560 | 139 | 0.837 | 0.248 | 0.383 | +0.000 | 0.569 | 0.566 | **0.457** | 0.132 | 0.338 |
| W3 E=250/500 | 560 | 139 | 0.837 | 0.248 | 0.383 | +0.000 | 0.569 | 0.566 | 0.410 | 0.132 | 0.338 |
| W3 E=1000/2000 | 560 | 139 | 0.837 | 0.248 | 0.383 | +0.000 | 0.569 | 0.566 | 0.419 | 0.132 | 0.338 |
| W4 | 557 | 139 | 0.837 | 0.250 | 0.385 | +0.002 | 0.569 | 0.566 | 0.404 | 0.132 | 0.338 |
| W5 | 557 | 139 | 0.837 | 0.250 | 0.385 | +0.002 | 0.569 | 0.566 | 0.459 | 0.132 | 0.338 |
| **ANN ceiling** | 166 | 166 | 1.000 | 1.000 | **1.000** | +0.617 | 1.000 | 1.000 | 1.000 | 1.000 | 0.839 |

Boundary offsets on (b), W0: median 5′ offset **+1875 bp** (the node starts 1875 bp *inside* the gene's 5′ end), median |5′| 2564; median 3′ offset −3 bp, median |3′| 505. **The de novo node is 5′-truncated, not 3′-truncated.** On (a) and (c) the medians are 0 and +150 / −1.

### (c) gorilla, held back, locus level only — 2664 nodes, 1975 expressed records

| variant | nodes | matched | R | P | **F** | ΔF | pair recall | pair prec |
|---|---|---|---|---|---|---|---|---|
| W0 | 2664 | 1824 | 0.924 | 0.685 | **0.786** | — | 0.683 | 0.906 |
| W1 / W2 / W3 (whole grid) | 2664 | 1824 | 0.924 | 0.685 | 0.786 | +0.000 | 0.683 | 0.906 |
| W4 (138 folds) | 2526 | 1818 | 0.921 | 0.720 | **0.808** | +0.021 | 0.684 | 0.907 |
| W5 | 2526 | 1818 | 0.921 | 0.720 | 0.808 | +0.021 | 0.684 | 0.907 |
| **ANN ceiling** | 1975 | 1975 | 1.000 | 1.000 | 1.000 | +0.214 | 1.000 | 1.000 |

## 3. Metric 3 — FAMILY level on the NPIP truth, substrate (b), frozen scorer, triangle leaders

minimap2 was re-run only for query md5 keys absent from the captured union PAFs: 2149 new tx md5 (5.2 Mbp, 110 s) and 1557 new body md5 (50.8 Mbp, 3 batches, 88+48+38 s), against the prebuilt `npip_ladder/idx/target.{splice,asm20}.mmi`.

| arm | nodes | exon / body edges | pairs | fams | **R** | P NPIP-only | F | P strict | F strict | full-length | mean copy cov |
|---|---|---|---|---|---|---|---|---|---|---|---|
| W0 | 560 | 1349 / 1885 | 1919 | 68 | **1.000** | 1.000 | 1.000 | 0.290 | 0.450 | 5/27 | 0.517 |
| W1 | 560 | 1349 / 1885 | 1919 | 68 | 1.000 | 1.000 | 1.000 | 0.290 | 0.450 | 5/27 | 0.517 |
| **W2** | 560 | 1349 / **4365** | 4761 | 22 | **0.778** | 1.000 | 0.875 | 0.273 | 0.404 | 5/27 | 0.517 |
| W3 E=500 | 560 | 1349 / 2068 | 2105 | 66 | 1.000 | 1.000 | 1.000 | 0.278 | 0.435 | 5/27 | 0.517 |
| W3 E=2000 | 560 | 1349 / 2543 | 2598 | 57 | 1.000 | 1.000 | 1.000 | 0.270 | 0.425 | 5/27 | 0.517 |
| **W4** | 557 | 1328 / 1864 | 1896 | 68 | 1.000 | 1.000 | 1.000 | **0.297** | **0.458** | 5/27 | 0.517 |
| W5 | 557 | 1328 / 4341 | 4728 | 22 | 0.778 | 1.000 | 0.875 | 0.276 | 0.408 | 5/27 | 0.517 |
| W7k3 (post hoc) | 560 | 2800 / 4127 | 5168 | 31 | 1.000 | 1.000 | 1.000 | 0.130 | 0.231 | **13/27** | 0.732 |
| **W9k3 (post hoc)** | 560 | 1796 / 2400 | 2485 | 65 | 0.889 | 1.000 | 0.941 | **0.480** | **0.623** | 12/27 | 0.673 |
| W10 (post hoc) | 557 | 1769 / 4486 | 5128 | 20 | 0.889 | 1.000 | 0.941 | 0.296 | 0.444 | 12/27 | 0.673 |
| **ANN ceiling** | 166 | 455 / 678 | 732 | 23 | 0.963 | 1.000 | 0.981 | 0.667 | 0.788 | 27/27 | 1.000 |

Two things to take from this table. **(i) W2 is not safe.** Putting terminal introns in the body more than doubles the gene-body edges (1885 → 4365), collapses 68 families into 22, and six of the 27 NPIP copies fall out of the matched family: FAMILY R 1.000 → 0.778. W5 inherits it. The same applies, more mildly, to the fixed extension: E = 2000 keeps R at 1.000 but costs P strict 0.290 → 0.270. **(ii) The only declared variant that improves anything at family level is W4**, which is also the only one that is free everywhere.

## 4. Post-hoc variants (labelled post hoc; not in the declarations)

Three were added after the declared curve came back flat, to find where the movement actually is.

- **W6 READ-BRIDGE MERGE** — nodes whose exons are both covered by ≥ 2 of the same reads are one locus (the ≥ 2-read rule already in `split_linked`; no new threshold). **Kill.** It over-merges: (b) R 0.837 → 0.446, 65 records lose their node; (a) R 0.960 → 0.920, 28 lost; (c) R 0.924 → 0.826, 192 lost. F falls on (a) and (b). Do not re-propose read-bridging as a locus-merge rule.
- **W7 READ-ISOFORM LOCUS** — the pieces of a node are the distinct intron chains of its assigned reads (≥ k reads each, blocks at the A3b median-low start/end); exons = merge(original ∪ chains); tx query set = those chains + the shipped rep, max over queries, never concatenated; body = first..last read base. This is W1 done from reads instead of from a piece list that is empty.
- **W9 JUNCTION-ANCHORED READ-ISOFORM LOCUS** — W7 but a chain is admitted only if it shares ≥ 1 junction with a chain already in the node (single-exon chains only if contained in the current exon span). It has to corroborate the node's splice structure. No distance threshold.

**On the ideal substrate the width fix works and is free:**

| | W0 | W7k3 | W9k3 | ANN |
|---|---|---|---|---|
| per-locus recall (a) | 0.858 | **0.993** | 0.987 | 1.000 |
| per-locus precision (a) | 0.991 | 0.990 | 0.990 | 1.000 |
| transcripts fully inside the node (a) | 0.603 | **0.975** | 0.960 | 1.000 |
| NPIP copies with a full-length node (a) | **11/26** | **23/26** | 22/26 | 26/26 |
| NPIP mean copy coverage (a) | 0.611 | 0.909 | 0.876 | 1.000 |
| TBC1D3 full-length (a) | 14/16 | 16/16 | 16/16 | 16/16 |

That is the established 11/27 figure, reproduced (11/26 here; one copy is off chr16/17/18 in the ideal locus set), and 12 of the 15 missing full-length copies recovered — at a per-locus precision cost of 0.001.

**On real testis reads it does not transfer.** (b): per-locus precision 0.566 → 0.342 (W7k2) / 0.373 (W7k3) / 0.436 (W9k3); junction presence actually *falls* 0.338 → 0.268 (W7k2), because unspliced pre-mRNA reads and 5′-run-on chains are admitted as "isoforms". Bipartite F rises only to 0.405. Gorilla, the held-back substrate, agrees with testis and not with the ideal arm: pair precision 0.906 → 0.828 (W7k3) / 0.871 (W9k3), F unchanged.

The junction anchor is what makes it survivable: on (b), W9k3 keeps FAMILY **P strict 0.480 and F strict 0.623** against W0's 0.290 / 0.450 and the annotated ceiling's 0.667 / 0.788 — the largest family-level gain measured in this experiment — with 12/27 full-length copies instead of 5/27. It costs FAMILY R 1.000 → 0.889 (3 copies land outside the matched family). W10 (W9 + W4 + W2) throws the gain away again through W2's body extension (P strict back to 0.296), which is a second, independent demonstration that W2 is the harmful ingredient.

## 5. Decision-rule outcome (verbatim rule in §4 of the declarations)

"rises by ≥ 0.05 ... on substrate (a) and (b)": largest ΔF on (a) is **+0.001** (W4), on (b) **+0.022** (post-hoc W7k3). **No variant, declared or post hoc, closes the gap.** Clause by clause:

| variant | ΔF (a) | ΔF (b) | Δ pair precision (b) | FAMILY R (b) | records losing their node (a/b) | verdict |
|---|---|---|---|---|---|---|
| W1 | +0.000 | +0.000 | +0.000 | 1.000 | 0 / 0 | fails (no gain; no-op) |
| W2 | +0.000 | +0.000 | +0.000 | **0.778** | 0 / 0 | fails (no gain, and R falls) |
| W3 (all E) | +0.000 | +0.000 | +0.000 | 1.000 | 0 / 0 | fails (no gain) |
| W4 | +0.001 | +0.002 | +0.000 | 1.000 | 0 / 0 | fails the ≥ 0.05 bar; passes every guard |
| W5 | +0.001 | +0.002 | +0.000 | **0.778** | 0 / 0 | fails (no gain, and R falls) |
| W6 post hoc | −0.018 | −0.029 | −0.216 | not scored | 28 / 65 | fails everything |
| W7k3 post hoc | −0.001 | +0.022 | −0.192 | 1.000 | 10 / 2 | fails |
| W9k3 post hoc | −0.001 | −0.006 | −0.129 | **0.889** | 7 / 2 | fails |

§6 of the declarations named this outcome in advance: **"locus width as parameterised here does not close the gap"** is the result, on the declared metric. But the declared metric is the wrong instrument for two of the three substrates, and that has to be said plainly:

- On (a) the bipartite F ceiling is 1.000 and W0 already sits at **0.954**. There are only 0.046 points of headroom, so a ≥ 0.05 rise is arithmetically impossible. The deficit on (a) lives entirely in the *per-pair* quantities (recall 0.858, isoform containment 0.603), and W7/W9 close those almost completely.
- On (b) the bipartite F is pinned by node count, not by node width. **261 of the 560 nodes overlap no annotated record's exons at all, and 383 of 560 are single-exon.** With 166 expressed records the precision ceiling on that node set is 166/560 = **0.296**, so F cannot exceed ≈ 0.44 no matter how the nodes are shaped. The 5′ truncation (median +1875 bp) and the 0.569 per-pair recall are genuine width defects, but they are second order behind node admission.

## 6. What to do with this

1. **Ship W4 (fragment fold) — it is free.** Never harmful anywhere: locus F +0.001/+0.002/+0.021 on (a)/(b)/(c), FAMILY P strict 0.290 → 0.297, F strict 0.450 → 0.458, R held at 1.000, 0 records lost on (a) and (b) (6 on gorilla, which needs a look before it goes in: gorilla folds 138 nodes where human folds 3).
2. **Do not ship W2 or W5, and do not ship any body extension.** The terminal-intron body is the one change that breaks the family (R 1.000 → 0.778); the fixed extension erodes P strict monotonically in E with no compensating gain.
3. **W1 as written cannot help, and the reason is upstream.** The node has one piece because the catalog picks one rep per locus. If the goal is "loci that encompass all isoforms of a gene", the isoforms have to be re-derived from the reads at the node (W7/W9) or the catalog has to stop collapsing to one rep — changing `shared_definition.rs` alone cannot do it.
4. **W9 (junction-anchored read-isoform locus) is the one worth developing.** On the ideal substrate it is nearly free and recovers 22/26 NPIP full-length copies; on real reads it is the only construction that raises FAMILY F strict (0.450 → 0.623) rather than lowering it. Its open problem is the FAMILY R cost (1.000 → 0.889), which is a family *split*, not an over-merge, and should be diagnosed per copy before anything is pre-registered.
5. **Do not re-propose read-bridge merging (W6).** Clean kill on all three substrates.

## 7. Provenance

- Substrate (a) `/mnt/linuxdisk/home/juanfraitu/npip_ideal/` (bam/ideal.bam, nodes/nodes.IDEAL.tsv, logs/rust2.err, transcripts.tsv, locus_set.tsv), chr16/17/18, 129 690 primary MAPQ ≥ 1 reads.
- Substrate (b) `/mnt/linuxdisk/home/juanfraitu/npip_ladder/` (reads.bam, rust/dump/x.nodes.tsv, union/*.paf, scripts/ladder.py for the mirror and the frozen scorer), 52 917 reads, 18 windows.
- Substrate (c) `/mnt/linuxdisk/home/juanfraitu/ggo_npip/` (dn/GGO.3ctg.bam, nodes/GGO.dn.nodes.tsv, nodes/GGO.3ctg.gff), 139 912 reads, 3 contigs.
- New alignments: `locus_width/map/tx_b00.paf`, `body_b0{0,1,2}.paf`, produced with the shipped flags `minimap2 -c -N 50 -p 0.1 -x splice -uf -t 4` and `-c -N 50 -p 0.1 -x asm20 -t 4` against the prebuilt `npip_ladder/idx/*.mmi`. Every other query key was served from the captured PAFs.
## Verification (independent recompute) -- agent 2 of 2

My code: `/mnt/linuxdisk/home/juanfraitu/locus_width/verify/` (`vr_gff.py` GFF/record builder, `vr_reads.py` BAM -> strand-corrected exon blocks, `vr_metrics.py` Jaccard + component-wise `linear_sum_assignment`, `vr_variants.py` W2/W3/W4, `vr_iso.py` metric 2, `vr_w7.py` my own read-isoform construction, `vr_all.py` -> `vr_SUMMARY.tsv`). I did not open any file under `locus_width/scripts/`; from the builder's workspace I read only `DECLARATIONS.txt`, `out/SUMMARY_locus.tsv`, `out/locus_{a,b,c}.json`, `out/family_b.json`, `out/sub_b.pkl`, `out/key2md5_b.pkl` and `map/*.paf`.

**1. Declarations precede every result file. CONFIRMED.** `DECLARATIONS.txt` mtime 2026-09-17 19:16:07 local, md5 `ec5a4a2812186e266e6395b212614e9b` (matches the report), `DECLARATIONS.time` = 2026-09-18T02:16:07Z. Earliest builder script 19:17:44, earliest output 19:20:22, `SUMMARY_locus.tsv` 19:36:11. The declarations fix W0-W5 + ANN, the four decision-rule clauses and their operational reading; W6-W10 are correctly labelled POST HOC per its own clause "any further variant is labelled POST HOC".

**2. Node sets rebuilt and diffed node by node.**
- (b): builder `sub_b.pkl['nodes']` is field-for-field identical (chrom, strand, n_reads, exons, rep_exons) to `npip_ladder/verify/nodes.pkl['arms']['A4']`, the node set produced by the *previous, independent* verification round -- 560/560, zero diffs. Its `body` equals the exon-union span for 560/560. Reads (52,917), regions (18) and the 27-copy NPIP truth are also bit-identical to that verified mirror. **Parity claim confirmed.**
- (a): 710 nodes on chr16/17/18 in `nodes.IDEAL.tsv`, exactly as claimed; `body` == exon-union span for 710/710.
- (c): 2664 nodes in `GGO.dn.nodes.tsv`; body == exon span for 2631 (33 nodes carry a wider table span -- harmless, 1b only).
- Variant node sets: W2/W3/W4/W5 rebuilt from scratch from my own read assignment. W4 reproduces 557 (b) and 2526 (c) exactly; (a) does not (see correction 1).

**3. Truth denominators rebuilt from the raw inputs. All three match exactly.** From `chm13v2.0_RefSeq_full.gff.gz` + my own `samtools view -F 2308 -q 1` extraction with the ts:A:- strand flip: (b) 239 records in the 18 windows, **166 expressed** (>= 3 strand-matched reads on exons) -- identical id set to the builder, 0 exon-union diffs, 0 read-count diffs. (a) from `locus_set.tsv`/`transcripts.tsv` on chr16/17/18: 759 records, **701 expressed**. (c) from `GGO.3ctg.gff` + `GGO.3ctg.bam`: 4477 genes with exons, **1975 expressed**. The denominators are identical across every variant row of SUMMARY_locus.tsv, i.e. they are not conditioned on any variant's output.

**4. Metric recompute. Every declared row on (b) and (c) reproduces to 3 dp; (a) reproduces except W4.**
- (a) W0: nodes 710 / recs 701 / match 673, R 0.960 P 0.948 F 0.954, pairR 0.858 pairP 0.991, nodes_ge2rec 61, unmatched nodes 37 -- all exact. W1/W2/W3 grid identical to W0 on 1a (as the declarations predicted); span-form W3 0.955/0.957/0.957/0.958/0.958 exact. ANN 701/701/701 = 1.000 everywhere, ge2 89 exact.
- (b) W0: 560/166/139, R 0.837 P 0.248 F 0.383, pairR 0.569 pairP 0.566, recs_no_node 27, recs_no_ov 19, ge2 42, spanF 0.402 -- all exact. W2 spanF 0.457 (spanR saturates at 1.000), W3_E250/500 0.410, E1000/E2000 0.419, W4 0.385 / span 0.404, W5 0.385 / span 0.459 -- all exact. ANN 166/166 = 1.000, ge2 47 exact.
- (c) W0: 2664/1975/1824, R 0.924 P 0.685 F 0.786, pairR 0.683 pairP 0.906, ge2 248, unmatched 840, spanF 0.790 -- exact. W4 2526 nodes, 0.921/0.720/0.808, ge2 239, **6 records lose their node** -- exact. W3 span 0.791/0.793/0.794 exact. ANN 1.000, ge2 339 exact.
- Metric 2 tx_contained: (a) 0.602 vs 0.603, (b) 0.132 exact, (c) 0.076 exact, ANN 1.000 on all three. junc_present does not reproduce (correction 3).
- Family level (b): `family_b.json['W0']` reproduces the frozen A4 row of `npip_ladder/results.json` field for field -- 1349 exon edges, 1885 body edges, 1919 pairs, 68 families, 364 loci, FAMILY R 1.000 / Pm 1.000 / Ps 0.2903 / Fs 0.4500, full_length 5, mean_copy_cov 0.5172. W1 is bit-identical to W0 in every family field (no-op confirmed; the "12 extra tx keys" is arithmetically consistent with 8 two-piece + 2 three-piece nodes). W2 body edges 1885 -> 4365, families 68 -> 22, FAMILY R 1.000 -> 0.778 (= 21/27) confirmed from the stored rows; W9k3 Ps 0.480 / Fs 0.623 / FL 12, W7k3 Ps 0.130 / FL 13 confirmed. I did not re-run minimap2, so the family rows are verified as stored + baseline parity, not re-derived from alignments.

**5. Decision rule, applied independently. The builder's verdict holds.** Clause 1 (locus F +>= 0.05 on both (a) and (b)): the largest declared-variant delta in my recompute is **+0.002** (W4 on (b)) and **0.000** on (a) -- every declared variant fails by a factor of ~25. Clause 2 (pairP drop <= 0.05): all declared variants move pairP by 0.000. Clause 3 (no record loses its node): declared variants lose 0 records on (a) and (b). Clause 4 (FAMILY R must not fall): W2 and W5 fall 1.000 -> 0.778 and are disqualified; W1/W3/W4 hold. **"No declared variant closes the gap" is confirmed, and W2 must not ship.**

**6. Anti-traps.**
- *Denominator conditioned on own output*: clean. R's denominator is the fixed expressed-record set (166/701/1975, identical across all rows, built from reads + annotation only). P's denominator is the node count, which is the prediction -- standard precision, but it means a variant that only DELETES nodes gains F for free; W4 is exactly that case (correction 2).
- *Ceiling arm scores 1.000 against itself*: confirmed independently on all three substrates for 1a, 1b and tx_contained. Note the 1a ceiling is self-scoring by construction (nodes ARE the records) so it carries no information about achievable precision; and on FAMILY R the ANN arm scores 0.963, below W0's 1.000, so "ceiling" is the wrong word there.
- *Merges penalised, not rewarded*: confirmed structurally and empirically. One-to-one matching leaves the swallowed record unmatched, so R falls -- visible in W6 on (b) (R 0.837 -> 0.446) and (c) (1824 -> 1632 matched, 192 records lose their node). The `nodes_ge2rec` column is real: I recomputed it (node overlapping >= 2 records by >= 50 exonic bp each) and it matches the builder exactly on W0 and W4 for all three substrates (61/42/248 and 61/42/239), and on (a) W7 (64).
- *Gain only from bigger nodes swallowing genes*: yes on (b) for W7 (correction 7: ge2 42 -> 102/112 alongside the +0.022); no on (a) for W7 (ge2 61 -> 64 while pairP holds at 0.990); and the W4 gains everywhere are node deletion, not node growth.

**7. Independent test of the post-hoc headline.** Since W7/W9 are undefined in the declarations, I wrote my own read-isoform locus (node exons := merged union of the blocks of assigned reads whose exact exon chain has >= k supporting reads; junction-anchored variant keeps only multi-exon chains whose every junction has >= k support). On (a) it lands on the builder's numbers almost exactly: per-locus recall 0.858 -> **0.992** (builder 0.993), tx_contained 0.602 -> **0.972** (builder 0.975), per-locus precision 0.991 -> **0.992** (builder 0.990), matched 673 -> 672, ge2 61 -> **64** (builder 64), records losing their node **10** (builder 10). On (b) it confirms the over-extension: pairP 0.566 -> 0.34..0.55, ge2 42 -> 85..101, locus F only +0.005..+0.011. **The headline post-hoc finding -- rebuilding the node's exon set from the reads' isoforms closes the isoform-containment deficit on the ideal substrate at no precision cost, and over-extends on real testis reads -- is independently corroborated.** My (c) behaviour differs from the builder's (my W7 is flat there, my W9 raises pairR to 0.75 while ge2 goes 248 -> 399), which is expected from a different rule and is why W6-W10 need a written definition before they are cited.

My full independent table: `/mnt/linuxdisk/home/juanfraitu/locus_width/verify/vr_SUMMARY.tsv`.