# Per-copy recovery of NPIP and TBC1D3, re-scored on the CAT/Liftoff v2.0 truth (human, 2026-10-01)

This re-runs the scoring of `docs/COPY_RECOVERY_TOOLS_2026-09-29.md` (pre-registration
`docs/PREREG_copy_recovery_tools_2026-09-29.md`, sha1 40b581b8) with a truth built from the T2T-CHM13 v2.0 CAT/Liftoff
annotation, under `docs/CAT_RERUN_PROTOCOL_2026-10-01.md` (R1-R6 and Amendment 1). The registered RefSeq results stand
as registered; both are reported below with the same scorer.

- **Scope.** Human A119b only (NPIP on chr16, TBC1D3 on chr17). Gorilla is not touched.
- **What is re-run.** Truth construction, gffcompare, read support, per-copy model scoring, the trace port, and the
  report. The tool GTFs are the registered ones, unchanged (they were produced without any annotation).
- **Nothing is re-tuned (R6).** Every definition, denominator and threshold is the frozen one. Where CAT forced a choice
  the frozen code never met, the choice and its alternative are both reported (§7).
- **Status.** One run, by the session that wrote this. Not independently recomputed.

## 1. Answer first

1. **Sanity check passed exactly.** The adapted chain, run on the original RefSeq truth, reproduces the registered human
   numbers: byte-identical gffcompare tmaps/refmaps and read tables, identical per-copy model JSON, identical summary on
   every registered key (§2).
2. **Human NPIP moves a lot under CAT; human TBC1D3 barely moves.** On copies with at least 2 exact-chain reads (E2),
   human NPIP is: ours 8 of 12, StringTie 8, FLAIR 3, isoseq 10 (RefSeq: ours 11 of 13, StringTie 10, FLAIR 6,
   isoseq 9). Human TBC1D3 on E2 is 9 of 10 for all four methods (RefSeq: 9 of 9 for ours, StringTie and FLAIR; 8 of 9
   for isoseq).
3. **Most of the NPIP change comes from the annotation, not from the models.** The copies in E2 change. Two RefSeq
   E2 copies have CAT models whose chains no read carries: NPIPA2 has 140 exact-chain reads under RefSeq and 0 under
   CAT; NPIPA5 has 19 and 0. LOC124907834 falls from E2 to E0, NPIPB5 from E2 to E1, and NPIPB3 has no CAT gene. Four
   copies enter E2 through CAT transcripts RefSeq lacks: NPIPA1, the PKD1P6-NPIPP1 row, NPIPA9 and NPIPB14P. Every
   method loses NPIPA2 and NPIPA5. isoseq gains at NPIPA1 and NPIPA9, two CAT genes that also contain read-through
   transcripts.
4. **Applied to the CAT numbers, the pre-registered P3 fails and F1's condition is met at human NPIP:** isoseq beats
   ours on E2 by 2 copies (10 against 8). Of isoseq's 4 E2 copies that ours lacks:
   - NPIPA1: one 1-read 3-exon model matching a 2-junction CAT chain.
   - NPIPA9: a 5-junction chain of a 21-junction R3 read-through gene.
   - NPIPB6: the same 6-junction chain as under RefSeq, where ours is `j` because the polish ISM step drops it
     (r1064, as registered).
   - NPIPB7: CAT's chain differs from RefSeq's by one junction. Our 17-read model carries RefSeq's chain; isoseq's 2-
     and 3-read models carry CAT's.

   Ours has 2 E2 copies isoseq lacks (LOC128966608, LOC124907808). StringTie ties ours (8 against 8).
5. **Two problems in the truths, flagged rather than fixed:**
   - **The PKD1P6-NPIPP1 row is a PKD1P6 copy in both truths, not an NPIP copy.** The registered RefSeq territory
     (chr16:15,126,650-15,141,806) lies inside RefSeq's own PKD1P6 record (15,126,099-15,159,720). Step 1's R1 image is
     the Liftoff PKD1P6 gene LOFF_G0001012. An NPIPA5 mRNA maps to 15,105,370-15,124,455, which is CAT's NPIPP1
     (CHM13_G0020725). Both ours and isoseq are COMPLETE at this row only under CAT, and that COMPLETE is a PKD1P6 model.
     Sensitivity arm S2 (§7) puts NPIPP1 in this row: E2 becomes ours 7 of 11 against isoseq 9 of 11.
   - **CAT genes carry many short partial transcripts** (2-4 exons). COMPLETE and E2 then count matches to chains of 1-3
     junctions. Three CAT NPIP COMPLETE calls rest on such a chain: isoseq at NPIPA1 (2 junctions), isoseq at NPIPB4
     (3), StringTie at NPIPB4 (2). Ours has none; its shortest matched chain has 5 junctions. §6.2 has the per-copy
     depths.

## 2. Sanity check (adapted chain, original RefSeq truth)

Run: `gc_cat.sh` / `reads.py` (frozen copy, unchanged) / `models_cat.py … refseq` / `report_cat.py` against
`copy_recovery_tools/ann/truth.hsa.{gtf,json}` and the families_gw RefSeq partner tables, output in `sanity_refseq/`.

| check | result |
|---|---|
| gffcompare tmap and refmap, 5 tools | byte-identical to `copy_recovery_tools/gc/hsa.*.{tmap,refmap}` |
| `reads.hsa.json`, `reads.hsa.tsv` | byte-identical to the registered files |
| `models.hsa.<tool>.json` (copies, disagreements), 5 tools | equal to the registered files |
| `summary.json` (human) | equal on every registered key |
| human part of `report.md` | identical apart from the one added column |
| headline | NPIP COMPLETE /26: ours 11, StringTie 12, isoseq 9, FLAIR 6; on E2 /13: ours 11, StringTie 10, FLAIR 6, isoseq 9; TBC1D3 ours 11/11, E2 9/9; FUSED excl. readthrough-defined: ours 11, StringTie 12, FLAIR 15, isoseq 18 (= the registered bracketed counts) |

The tool GTFs are byte copies of the registered restricted files (sha1 302ccb37 ours, 5c95ab87 ours primaries-only,
bdd7fe61 StringTie, fa29fa73 FLAIR, 4e2f3e91 isoseq).

## 3. The truth under CAT: what changed

| item | RefSeq (registered) | CAT (this run) |
|---|---|---|
| how a copy is chosen | Dishuck table rows by GeneID (NPIP); RefSeq `description` "TBC1 domain family member 3" (TBC1D3) | the `cat_gene_id` of each row in `lit_subclusters_npip_dishuck_check.CAT.tsv` / `tbc1d3_members.CAT.tsv` (R1/R2; never a name) |
| NPIP copies | 26 chr16 (+ NPIPB1P on chr18, not scored) | **25** chr16 (+ NPIPB1P → CHM13_G0025964, not scored). RefSeq NPIPB3 (chr16:21,337,400-21,360,419) is dropped: CAT has no gene on that span |
| TBC1D3 copies | 16 chr17 records | 16 CAT genes (all re-keyed) |
| isoforms | the record's RefSeq transcripts; NPIPB14P borrowed PDXDC2P-NPIPB14P's exons; PKD1P6-NPIPP1 clipped to exons ≥ 15,126,650 | all CAT transcripts of the gene. Both special cases are gone: NPIPB14P → CHM13_G0022125 has 3 transcripts of its own, and LOFF_G0001012 is used whole |
| scorable for COMPLETE | NPIP 26/26; TBC1D3 11/16 (5 records without exon features) | NPIP **24/25**; TBC1D3 **15/16** (§3.1) |
| read-through flag | 2 readthrough-defined copies (PKD1P6-NPIPP1, NPIPB14P) | 3 copies by R3, each as the R1 image of a RefSeq read-through: PKD1P6-NPIPP1 → LOFF_G0001012, NPIPA6 → LOFF_G0001021 (image of LOC131696449 = PKD1P1-NPIPA5L), NPIPA9 → CHM13_G0020804 (image of PKD1P5-LOC105376752). No truth gene has an A-B name |
| partner records (FUSED) | families_gw RefSeq `genes.tsv` + `genes_only.gff` (incl. RefSeq read-throughs) | `chm13v2.0_CAT_Liftoff.genes.tsv`, every CAT gene on chr16/chr17, incl. read-throughs |
| truth transcripts | 147 human chains | 179 chains (NPIP 25 copies, TBC1D3 16), 0 single-exon |

### 3.1 Single-exon and all-non-canonical CAT transcripts

The frozen construction merges every non-canonical intron. It had no single-exon chain at all ("0 single-exon truth
chains"). Under CAT, 13 transcripts at 12 copies end up with no junction:

- 10 are genuinely single-exon (or 2-4-exon with every intron non-canonical) and sit beside multi-exon transcripts of
  the same gene.
- **NPIPB13 → LOFF_G0001077** (Liftoff, named NPIPA3) has 3 transcripts whose 18 introns are all non-canonical
  (GT-TT, GA-TT, AT-TC, …). No constant coordinate shift repairs them; the best shift (−2 bp) makes 3 of 7 canonical.
- **TBC1D3P7 → CHM13_G0024278** (8 exons, 7 introns: CA-AC, AG-CA, TT-GG, …) is the same case.

These transcripts stay in the territory but are **not entered as chains**, so NPIPB13 and TBC1D3P7 are not scorable
(`n/s`). If they were entered, as the frozen code would do, every unspliced own read would count as an exact-chain
read and every single-exon model would be a candidate `=`. Arm S1 (§7) shows the effect.

Four of the five RefSeq TBC1D3 records without exon features are now scorable: TBC1D3P4, TBC1D3P3, TBC1D3P1, and
LOC100420311. LOC100420311's image is CAT TBC1D29P (CHM13_G0023800, span Jaccard 0.258).

### 3.2 Territories and chains change a great deal

CAT's gene models differ from RefSeq's at most NPIP copies:

- **Truncated models.** NPIPA2 has 2 transcripts with 4 junctions (RefSeq: 14 transcripts, 11-junction primary), and
  its territory drops from 4,517 to 820 bp. NPIPA5 has 3 transcripts (RefSeq: 12).
- **Read-through-length transcripts inside the gene.** The NPIPA1 gene has a 27-exon transcript, NPIPA6 25, NPIPA9 23,
  LOC128966608 23, NPIPB8 → AC138894.1 up to 23, and LOFF_G0001012 24. Territories grow accordingly; NPIPA1 goes from
  1,084 to 7,342 bp.
- **Many 2-4-exon partial transcripts** (§6.2).
- **The three NPIPB15 copies** each have a single multi-exon chain (RefSeq: 6-7).

CAT names are permuted against RefSeq. For example, RefSeq NPIPB5 is CAT "NPIPB3", RefSeq NPIPB13 is Liftoff "NPIPA3",
RefSeq TBC1D3F is Liftoff "TBC1D3B", and RefSeq TBC1D3E/K/D are CAT K/D/E. All joins here go by cid and gene id.

### 3.3 R2: name-based CAT sets reported beside the R1 image (not substituted)

**NPIP.** CAT genes on chr16 whose `gene_name` starts with NPIP and has no `-`: 23 genes.

- 21 of them are in the R1 truth.
- 2 are not: **NPIPP1** (CHM13_G0020725, 15,105,363-15,124,458, −) and **NPIPB7** (CHM13_G0021068,
  28,924,554-28,939,395, +; nested in NPIPB8's image).
- The R1 truth has 4 genes without an NPIP name: LOFF_G0001012 (PKD1P6), LOFF_G0001021 and CHM13_G0020804 (both
  AC138969.1), and CHM13_G0021067 (AC138894.1).

**TBC1D3.** CAT genes on chr17 whose name starts with TBC1D3:

- They are the 15 named images plus **TBC1D3J** (CHM13_G0024042, nested in TBC1D3B's span).
- The R1 set adds TBC1D29P, the image of LOC100420311.
- Elsewhere in the genome the prefix also matches TBC1D3P6 (chr1) and TBC1D30/31/32. Those are different genes, listed
  only because the rule matches the prefix.

The full list is in `ann/r2_namebased.hsa.tsv`.

## 4. Headline, side by side (RefSeq registered → CAT)

Human A119b — **NPIP**. RefSeq: 26 copies, 26 scorable, E2 13, E1 2. CAT: 25 copies, 24 scorable, E2 12, E1 2. Every
copy has at least 2 own primaries in both truths, so the "≥ 2 primaries" columns equal the "all" columns.

| arm | LOCUS | COMPLETE (all scorable) | COMPLETE on E2 | on E1 | partial-or-better (=,c,k) | FUSED copies | FUSED excl. read-through-flagged | UNIQUE | models per LOCUS copy, median / max |
|---|---|---|---|---|---|---|---|---|---|
| ours (default) | 26/26 → 25/25 | 11/26 → **9/24** | 11/13 → **8/12** | 0/2 → 1/2 | 15/26 → 11/24 | 13 → 8 | 11 of 24 → 6 of 22 | 13 → 17 | 16/76 → 14/76 |
| ours, primaries only | 26/26 → 25/25 | 10/26 → 8/24 | 10/13 → 8/12 | 0/2 → 0/2 | 12/26 → 11/24 | 12 → 6 | 10 of 24 → 5 of 22 | 14 → 19 | 7/42 → 8/42 |
| StringTie 3.0.1 -L | 26/26 → 25/25 | 12/26 → **8/24** | 10/13 → **8/12** | 2/2 → 0/2 | 15/26 → 10/24 | 14 → 12 | 12 of 24 → 10 of 22 | 12 → 13 | 13/37 → 13/37 |
| FLAIR 3.0.1 | 26/26 → 25/25 | 6/26 → **3/24** | 6/13 → **3/12** | 0/2 → 0/2 | 12/26 → 10/24 | 17 → 14 | 15 of 24 → 11 of 22 | 9 → 11 | 20/92 → 16/93 |
| isoseq collapse | 26/26 → 25/25 | 9/26 → **10/24** | 9/13 → **10/12** | 0/2 → 0/2 | 21/26 → 20/24 | 20 → 18 | 18 of 24 → 15 of 22 | 6 → 7 | 74/321 → 75/325 |

Human A119b — **TBC1D3**, 16 copies in both truths. RefSeq: 11 scorable, E2 9, E1 0. CAT: 15 scorable, E2 10, E1 0.
No copy is read-through-flagged in either truth.

| arm | LOCUS | COMPLETE (all scorable) | COMPLETE on E2 | partial-or-better | FUSED | UNIQUE | models median / max |
|---|---|---|---|---|---|---|---|
| ours (default) | 13/16 → 13/16 | 11/11 → **11/15** | 9/9 → **9/10** | 11/11 → 12/15 | 5 → 5 | 8 → 8 | 12/21 → 12/21 |
| ours, primaries only | 13/16 → 13/16 | 9/11 → 9/15 | 9/9 → 9/10 | 9/11 → 10/15 | 6 → 5 | 7 → 8 | 6/20 → 6/21 |
| StringTie | 14/16 → 14/16 | 9/11 → 9/15 | 9/9 → 9/10 | 9/11 → 10/15 | 5 → 4 | 9 → 10 | 12/41 → 11/43 |
| FLAIR | 13/16 → 13/16 | 9/11 → 9/15 | 9/9 → 9/10 | 9/11 → 11/15 | 9 → 8 | 4 → 5 | 14/39 → 14/40 |
| isoseq | 16/16 → 16/16 | 8/11 → 9/15 | 8/9 → 9/10 | 9/11 → 11/15 | 8 → 6 | 8 → 10 | 44/166 → 44/191 |

gffcompare `=` and the Python exact-chain check agree on every model in every arm (0 disagreements). No arm has a model
at a copy with 0 own primaries, because there is no such copy in either truth.

**First losing step (ours, absent copies).** The trace port was re-run on the CAT territories. The same three TBC1D3
pseudogenes are absent, with the same steps as registered:

- TBC1D3P4: "polish (gate survivors exist)"; the registered no-polish arm named this the mono shadow.
- TBC1D3P3: gate (strict junctions + single-exon `+` placeholder).
- TBC1D3P7: "polish"; registered as the mono floor.

The trace counts (own primaries 2/4/2, gate survivors 1/0/1) are unchanged. The no-polish arm was not re-run.

## 5. Per copy: what changed and why

Cell = best class against the copy (`=` is COMPLETE), shown as RefSeq → CAT when they differ. The exact-chain column
gives own primaries carrying an annotated chain of the copy, with the read class.

**NPIP**

| cid | RefSeq row | CAT gene (CAT name; R1 quality) | exact-chain reads R / C | ours | ours, primaries | StringTie | FLAIR | isoseq |
|---|---|---|---|---|---|---|---|---|
| h00 | NPIPB2 | CHM13_G0020650 (NPIPB2; strong) | 7 (E2) / 29 (E2) | j → **=** | j → **=** | j → **=** | **=** | **=** |
| h01 | NPIPA2 | CHM13_G0020702 (NPIPA2; partial) | 140 (E2) / 0 (E0) | **=** → j | **=** → j | **=** → j | **=** → j | **=** → j |
| h02 | NPIPA1 | CHM13_G0020714 (NPIPA1; partial) | 0 (E0) / 10 (E2) | j → c | j → c | j → c | j → c | c → **=** |
| h03 | PKD1P6-NPIPP1 | LOFF_G0001012 (PKD1P6; strong) | 0 (E0) / 4 (E2) | j → **=** | j → **=** | j | j | j → **=** |
| h04 | NPIPA5 | CHM13_G0020732 (NPIPA5; partial) | 19 (E2) / 0 (E0) | **=** → j | **=** → j | **=** → j | **=** → j | **=** → c |
| h05 | NPIPA6 | LOFF_G0001021 (AC138969.1; partial) | 0 (E0) / 1 (E1) | k → j | j | k → j | k → j | c |
| h06 | NPIPA7 | CHM13_G0020773 (NPIPA7; strong) | 0 / 0 | j | j | j | j | c |
| h07 | NPIPA8 | LOFF_G0001030 (NPIPA8; strong) | 0 / 0 | j | j | j | j | j → c |
| h08 | NPIPA9 | CHM13_G0020804 (AC138969.1; partial) | 0 (E0) / 11 (E2) | k → j | k → j | j | k → c | c → **=** |
| h09 | NPIPB3 | none: dropped | 8 (E2) / - | **=** → - | **=** → - | k → - | j → - | c → - |
| h10 | LOC128966608 | CHM13_G0020898 (NPIPB5; partial) | 8 (E2) / 36 (E2) | **=** | **=** | **=** | j → c | c |
| h11 | NPIPB4 | CHM13_G0020921 (NPIPB4; partial) | 2 (E2) / 26 (E2) | **=** | **=** | **=** | j → c | **=** |
| h12 | NPIPB5 | CHM13_G0020937 (NPIPB3; weak) | 2 (E2) / 1 (E1) | **=** | c | **=** → j | c | c |
| h13 | NPIPB6 | CHM13_G0021048 (NPIPB6; partial) | 2 (E2) / 2 (E2) | j | j | **=** | c → j | **=** |
| h14 | NPIPB7 | CHM13_G0021054 (NPIPB8; partial) | 7 (E2) / 4 (E2) | **=** → j | **=** → j | **=** → j | **=** → j | **=** |
| h15 | NPIPB8 | CHM13_G0021067 (AC138894.1; partial) | 0 / 0 | j → c | j → c | j → c | j → c | j → c |
| h16 | NPIPB9 | CHM13_G0021074 (NPIPB9; strong) | 0 / 0 | c → j | j | c → j | j | c → j |
| h17 | NPIPB10P | CHM13_G0021099 (NPIPB10P; strong) | 0 / 0 | j | j | j | j | j |
| h18 | NPIPB11 | CHM13_G0021115 (NPIPB11; strong) | 0 / 0 | j | j | j | j | c |
| h19 | NPIPB12 | CHM13_G0021130 (NPIPB12; strong) | 1 (E1) / 0 (E0) | c → j | j | **=** → j | j | c → j |
| h20 | LOC124907834 | CHM13_G0021187 (NPIPB13; weak) | 3 (E2) / 0 (E0) | **=** → j | **=** → j | **=** → j | c | c |
| h21 | NPIPB13 | LOFF_G0001077 (NPIPA3; partial) | 1 (E1) / n/s | j → n/s | j → n/s | **=** → n/s | j → n/s | c → n/s |
| h22 | NPIPB14P | CHM13_G0022125 (NPIPB14P; span) | 0 (E0) / 11 (E2) | j → **=** | j → **=** | j → **=** | j | j → **=** |
| h23 | NPIPB15 | CHM13_G0022247 (NPIPB15; partial) | 3 (E2) / 70 (E2) | **=** | **=** | j → **=** | **=** | **=** |
| h24 | LOC124907808 | LOFF_G0001213 (NPIPB15; partial) | 4 (E2) / 2 (E2) | **=** | **=** | **=** | c → j | **=** → k |
| h25 | LOC124907807 | LOFF_G0001218 (NPIPB15; partial) | 12 (E2) / 11 (E2) | **=** | **=** | **=** | **=** | **=** |

What the rows say:

- **Lost by everyone with the annotation.**
  - NPIPA2: CAT's 4-junction chain shares 3 junctions with RefSeq's primary. 140 reads carry RefSeq chains and none
    carries CAT's.
  - NPIPA5: two of CAT's three chains are also RefSeq chains, but not the ones the 19 exact-chain reads carry (a
    7-junction chain with 17 reads and two 8-junction chains with 1 each).
  - NPIPB3: no CAT gene.
- **Lost by ours (and StringTie, FLAIR) through a one-junction difference.**
  - NPIPB7: CAT's chain uses intron 28,747,569-28,757,212 where RefSeq has 28,747,569-28,751,033. Our 17-read model
    carries RefSeq's chain; isoseq's 2-3-read models carry CAT's.
  - LOC124907834: CAT's chains share 3 of RefSeq's 6 junctions, and no read is exact under CAT.
- **Gained by ours through CAT transcripts.**
  - NPIPB2: our 9-read 6-exon model matches a 5-junction CAT chain. Under RefSeq, ours lost NPIPB2's chains to the
    fraction rule (r1063).
  - LOFF_G0001012: a 7-junction PKD1P6 chain (§1.5).
  - NPIPB14P: a 5-junction chain of the copy's own CAT gene. RefSeq NPIPB14P's only isoform was a 25-junction
    read-through.
- **Gained by isoseq**: NPIPA1 (2-junction chain, 1-read model), NPIPA9 (5-junction chain of a 21-junction read-through
  gene), LOFF_G0001012, NPIPB14P.
- **The E1 cell turns over.** The RefSeq E1 copies were NPIPB12 and NPIPB13, where StringTie's single-read models were
  `=`. Under CAT, the 7-junction RefSeq chain carried by NPIPB12's single read is not among CAT's 8 chains (0 exact
  reads), and NPIPB13 is not scorable. The CAT E1 copies are NPIPA6 and NPIPB5. Ours is `=` at NPIPB5, the same 6-junction model as under RefSeq;
  NPIPB5 drops from E2 to E1 because RefSeq's 7-junction chain is not in CAT.
- **FUSED.** Ours goes 13 → 8 because:
  - NPIPA1 and NPIPA9 are no longer FUSED: their CAT genes include the PKD1P3 / PKD1P5 read-through exons, so those
    exons now lie inside the territory.
  - NPIPB9 (EIF3C) and NPIPB13 (SMG1P5) lose their partners under the CAT territories.
  - NPIPB3 is dropped.

  StringTie becomes FUSED at NPIPB4 (SMG1P4) and LOC124907834, and stops being FUSED at NPIPA6, NPIPB8 and NPIPB9.

**TBC1D3.** All 11 originally scorable copies keep the same COMPLETE call for every method. One best class changes
(StringTie at TBC1D3B, n → m). Two things change:

- **The four newly scorable records.** TBC1D3P4, TBC1D3P3, TBC1D3P1 and the TBC1D29P image are `=` for no method,
  except isoseq at TBC1D3P1: a 12-junction full chain, 4 exact-chain reads, model support 2. That copy enters E2, so
  isoseq's E2 rises from 8/9 to 9/10 and the other methods' from 9/9 to 9/10.
- **FUSED.**
  - TBC1D3G loses its partner LOC101060212 for StringTie, FLAIR and isoseq.
  - TBC1D3P1 becomes FUSED for ours and StringTie, with the CAT read-through TBC1D3P1-DHX40P1; FLAIR and isoseq were
    already FUSED there.
  - LOC100420311 is no longer FUSED with TBC1D29P (ours, StringTie, isoseq), because TBC1D29P is now the copy itself.

The full per-copy table (territory bp, chain counts, own primaries, FUSED per tool) is `score/compare_refseq_cat.md`.

## 6. Interpreting the NPIP change

### 6.1 The pre-registered statements, re-read on CAT (descriptive; nothing re-tuned)

| registered statement | RefSeq | CAT |
|---|---|---|
| P3: on E2, ours ≥ every baseline in each testable cell | holds (11/13 vs 10, 6, 9; TBC1D3 9/9 vs 9, 9, 8) | **fails at human NPIP** (8/12 vs isoseq 10); holds at TBC1D3 (9/10 for all four) |
| F1: a baseline exceeds ours on E2 by ≥ 2 copies in a cell | not met | **met at human NPIP** (isoseq +2); not met at TBC1D3 |
| "StringTie's one-read floor recovers the single-read copies" (E1) | StringTie 2/2 (NPIPB12, NPIPB13) | StringTie 0/2; ours 1/2 (NPIPB5) |
| seeding: ours COMPLETE at TBC1D3B and TBC1D3I with 0 exact-chain primaries, primaries-only arm loses them | holds | holds (unchanged) |
| all methods FUSED at the read-through copies | 4/4 for every method | the 3 R3-flagged copies: ours FUSED at 2 (LOFF_G0001012, NPIPA6), StringTie 2, FLAIR 3, isoseq 3 |

### 6.2 How deep are the COMPLETE matches? (descriptive)

For each COMPLETE call, this compares the longest truth chain the tool's `=` models match against the copy's longest
chain. Script: `code/complete_depth_cat.py`.

| copy (CAT) | longest chain | ours | StringTie | FLAIR | isoseq |
|---|---|---|---|---|---|
| NPIPB2 | 9 | 5 | 5 | 9 | 9 |
| NPIPA1 | 26 (read-through transcript) | - | - | - | **2** |
| LOFF_G0001012 (PKD1P6) | 23 | 7 | - | - | 7 |
| NPIPA9 | 21 (read-through transcript) | - | - | - | 5 |
| LOC128966608 | 22 | 8 | 8 | - | - |
| NPIPB4 | 8 | 6 | **2** | - | **3** |
| NPIPB5 | 8 | 6 | - | - | - |
| NPIPB6 | 6 | - | 6 | - | 6 |
| NPIPB7 | 6 | - | - | - | 6 |
| NPIPB14P | 6 | 5 | 5 | - | 5 |
| NPIPB15 (3 copies) | 6 | 6, 6, 6 | 6, 6, 6 | 6, -, 6 | 6, -, 6 |

Under RefSeq every NPIP COMPLETE call matched a chain of at least 6 junctions. Under CAT, three rest on 2-3-junction
chains (bold). The E2 class also admits copies through short chains. At NPIPA1, 8 of the 10 exact-chain reads carry a
3-junction chain, 1 a 2-junction chain and 1 a 5-junction chain. At NPIPB4, 14 reads carry a 2-junction chain and 7 a
1-junction chain; NPIPB4 would still be in E2 on its 6-junction chain alone (2 reads), as under RefSeq.

**Descriptive count only, not a re-scored measure:** dropping the three short-chain calls would give COMPLETE on E2 of
ours 8, StringTie 7, FLAIR 3, isoseq 8.

### 6.3 What it means

**The per-copy comparison of methods is sensitive to which annotation defines the copies and their chains, at NPIP
far more than at TBC1D3.** TBC1D3's CAT models for the protein-coding and transcribed copies are 12-13-junction
full-length models close to RefSeq's (R1 quality: 9 strong, 2 partial, 5 span-mapped pseudogene records), and every
method's E2 count is within one copy of RefSeq's. NPIP's CAT models differ in three ways: they
are often partial (R1 quality "partial" or "weak" at 16 of 25 copies), they mix read-through-length and 2-4-exon
partial transcripts in one gene, and one row is not an NPIP copy at all.

The registered conclusion, "where the reads carry the whole chain, ours is COMPLETE at least as often as every
baseline", holds on the RefSeq truth and on the CAT TBC1D3 truth. It does not hold on the CAT human NPIP truth as
built under the protocol. There, isoseq leads by 2 copies net. Of the 4 copies isoseq has and ours lacks:

- 2 rest on a short chain (NPIPA1, 2 junctions) or on a partial chain of a read-through gene (NPIPA9).
- 2 are full-length chains: NPIPB6, lost by ours at ISM as under RefSeq, and NPIPB7, whose CAT chain differs from
  RefSeq's by one junction.

## 7. Sensitivity arms (reported beside the headline, never substituted)

| arm | NPIP scorable / E2 / E1 | NPIP COMPLETE, all scorable: ours / ours-prim / StringTie / FLAIR / isoseq | NPIP COMPLETE on E2 | TBC1D3 scorable / E2 | TBC1D3 COMPLETE on E2 (all four methods) |
|---|---|---|---|---|---|
| CAT headline | 24 / 12 / 2 | 9 / 8 / 8 / 3 / 10 | 8 / 8 / 8 / 3 / 10 | 15 / 10 | 9/10 |
| S1: single-exon CAT transcripts entered as chains (the frozen code applied literally) | 25 / 16 / 0 | 9 / 8 / 8 / 3 / 10 | 9 / 8 / 8 / 3 / 10 (of 16) | 16 / 11 | 9/11 |
| S2: PKD1P6-NPIPP1 row = CAT NPIPP1 (CHM13_G0020725) instead of the R1 image | 24 / 11 / 2 | 8 / 7 / 8 / 3 / 9 | 7 / 7 / 8 / 3 / 9 (of 11) | 15 / 10 | 9/10 |

- **S1 changes no COMPLETE count, only the read classes.** Unspliced reads become "exact-chain" reads, for example
  NPIPB4 26 → 388, LOC128966608 36 → 305 and NPIPB5 1 → 62. E2 grows by 4 NPIP copies and 1 TBC1D3 copy (TBC1D3P7, 2
  unspliced reads). This is why the headline does not enter single-exon chains.
- **S2.** The real NPIPP1 copy has 387 own primaries, 0 exact-chain reads, and the largest same-chain group is 42.
  Every method has a FUSED locus there (with PKD1P6 and the CAT PKD1P6-NPIPP1 read-through), and the best class is
  `j` for every method. Ours and isoseq both lose the PKD1P6 COMPLETE call. On E2, isoseq still leads ours by 2
  (9 against 7).

## 8. Deviations

1. **Human only**, as instructed; the gorilla truths and `report.py`'s gorilla branch are not run.
2. **`gc_cat.sh` does not re-filter the tool GTFs.** It scores byte copies of the registered restricted files (sha1s in
   §2) and checks that they hold only chr16/chr17. Re-filtering a file onto itself would truncate it.
3. **Single-exon (after canonicalisation) CAT transcripts are kept in the territory but not entered as chains** (13
   transcripts at 12 copies). The frozen code never met this case. The literal behaviour is arm S1.
4. **cids keep the frozen numbering** (NPIP rows in file order, TBC1D3 by RefSeq start), so one cid names the same
   RefSeq source row in both truths. h09 (NPIPB3) is unused under CAT.
5. **R3 needs RefSeq descriptions.** They are read from the frozen chr16+chr17 slice
   (`copy_recovery_tools/ann/hsa_chr16_17.gff`, read only). An R1 image lies on the same chromosome, so the slice covers
   every truth gene.
6. **Partner names** are `gene_name|gene_id`, because CAT names repeat. The frozen conversion of partner exon blocks,
   `(a − 1, b)` on a 0-based table, is kept in both arms so that the sanity arm stays identical. It widens each partner
   exon by 1 bp at its start, under RefSeq as under CAT.
7. **`report_cat.py` adds one column**: FUSED copies excluding read-through-flagged copies. On the RefSeq arm it
   reproduces the registered post-hoc counts (11/10/12/15/18). `posthoc.py` was not re-run: its qualifiers hard-code
   RefSeq copy names.
8. **The trace port was re-run on the CAT territories** (heavy lock, 94 s, 2.3 GB). The regional no-polish arm
   (`arm_hsa.sh`) was not re-run; the polish-step names for TBC1D3P4/P7 are the registered ones.
9. **Arms S1 and S2 and the descriptive §6.2 are additions** that the pre-registration does not contain. None changes a
   headline cell.
10. **Two step-1 issues are flagged, not corrected here:** the PKD1P6-NPIPP1 row's R1 image (§1.5, S2), and the
    LOC100420311 → TBC1D29P image (span Jaccard 0.258).

## 9. Exact commands

All in `/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat` (`N`). `L` = `bash
/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/tools/rlock.sh light`; `O` = the original dir (read only). Everything ran
in the foreground.

```bash
# setup: code copies (originals read-only in code/frozen_orig), registered restricted tool GTFs, AS table link
cp -p $O/code/* $N/code/frozen_orig/; for T in ours ours_primary stringtie flair isoseq; do cp -p $O/runs/hsa.$T.gtf $N/runs/; done
ln -s /mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv $N/runs/hsa.molecules.tsv
# sanity arm (original RefSeq truth)
$L bash -c "for T in ours ours_primary stringtie flair isoseq; do bash $N/code/gc_cat.sh hsa \$T $O/ann/truth.hsa.gtf $N/sanity_refseq/gc; done"
$L python3 $N/code/frozen_orig/reads.py hsa $O/ann/truth.hsa.json $N/sanity_refseq/score/reads.hsa
$L bash -c "for T in ...; do python3 $N/code/models_cat.py hsa \$T $N/runs/hsa.\$T.gtf $N/sanity_refseq/gc/hsa.\$T.tmap $N/sanity_refseq/score/models.hsa.\$T $O/ann/truth.hsa.json refseq; done"
$L python3 $N/code/report_cat.py $O/ann $N/sanity_refseq/score $O/score/trace.hsa.json
# CAT arm
$L python3 $N/code/build_truth_cat.py --out $N/ann
$L bash -c "for T in ...; do bash $N/code/gc_cat.sh hsa \$T $N/ann/truth.hsa.gtf $N/gc; done"
$L python3 $N/code/frozen_orig/reads.py hsa $N/ann/truth.hsa.json $N/score/reads.hsa                 # 20 s, 74 MB
$L bash -c "for T in ...; do python3 $N/code/models_cat.py hsa \$T $N/runs/hsa.\$T.gtf $N/gc/hsa.\$T.tmap $N/score/models.hsa.\$T $N/ann/truth.hsa.json cat; done"
RLOCK_WAIT=480 bash .../rlock.sh heavy python3 $N/code/trace_port_cat.py hsa $N/ann/truth.hsa.json $N/score/trace.hsa.json
$L python3 $N/code/report_cat.py $N/ann $N/score $N/score/trace.hsa.json
python3 $N/code/compare_cat.py $O/ann $N/sanity_refseq/score $N/ann $N/score $N/score/compare_refseq_cat.md
python3 $N/code/complete_depth_cat.py $N/ann $N/score
# sensitivity arms: same chain with these truths, outputs in sens_keepmono/ and sens_a4_npipp1/
$L python3 $N/code/build_truth_cat.py --out $N/sens_keepmono/ann --keep-mono
$L python3 $N/code/build_truth_cat.py --out $N/sens_a4_npipp1/ann --override GeneID:105369154=CHM13_G0020725
```

## 10. Files

`/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/` (scratch, no backup; `RUN.sha1` records the code and
output sha1s):

- `code/`:
  - `build_truth_cat.py` (dcba5ef6), `gc_cat.sh` (a8e86c92), `models_cat.py` (688fb50d), `report_cat.py` (9e2ff0cf),
    `trace_port_cat.py` (8724555f), `compare_cat.py`, `complete_depth_cat.py`.
  - `frozen_orig/`: read-only copies of the frozen instruments, sha1s as in `FROZEN.orig.sha1`. `reads.py` is run from
    here unchanged.
- `ann/`:
  - `truth.hsa.{json,gtf}` and `copies.hsa.tsv` (the CAT truth).
  - `dropped.hsa.tsv` (NPIPB3).
  - `r2_namebased.hsa.tsv`.
- `gc/`, `score/`:
  - the CAT arm: `report.md`, `summary.json`, `reads.hsa.*`, `models.hsa.*`, `trace.hsa.*`.
  - `compare_refseq_cat.md`.
- `sanity_refseq/`: the RefSeq reproduction.
- `sens_keepmono/`, `sens_a4_npipp1/`: the two sensitivity arms.

The original `copy_recovery_tools/` dir was not modified: its frozen sha1s verify, and no file in it is newer than
2026-09-30.
