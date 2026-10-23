# B. BAM-tag and alignment-feature fingerprints of a missing copy (RNA reads on CHM13 only)

**Angle:** what a missing copy leaves in AS / NM / de / ms / MD / CIGAR / clips / SA / MAPQ after O2-style assignment. **Reads used:** BRIEF;
prereg Results; EXPLORE §1-5; expected failures §0-F3; register §6 and rows 1329, 1330, 1347, 1348, 1350, 1357, 1362; memory O3 notes; vault notes
(DuploMap, Parascopy, Paraphase 2025, SDA, IsoCon, longcallR, Clair3-RNA, Deschamps, Pal, Guitart, Dishuck, Greshnova, Tomaszkiewicz, Yoo).
**No truth table was opened.** The §1 numbers come from the frozen `reads.tsv`, `columns.tsv` and `calls.tsv`, plus the run's `teach/*.bam`. ⚠ I
looked at these distributions, including families the BRIEF names, so **these two donors are now seen data for every fingerprint below.**
A verdict needs a fresh donor (EXPLORE §5.4.3).

## 0. Bottom line
1. **The three classes need three different signals.** No single tag separates them.
   - (i) An absorbed donor-only copy adds sequence that no reference copy explains. It is detectable.
   - (ii) A reference copy taught by a paralog's haplotype adds nothing (EXPLORE §5.2: PRAMEF7.1 looks like a present homozygous copy). Only a
     positive test can expose it: does any read carry evidence that only that copy can produce?
   - (iii) Normal reads carry ≤ 2 haplotypes per copy, and their leftover edits are sequencing error.
2. **Raw NM and de are the wrong currency.** In HG002 tied reads, the median of min-over-copies `indels` is 1 base and the median of `subs` is 0.
   NM and de are therefore dominated by HiFi indel error. Use residual substitutions and recurrent indel *events*, not NM or de.
3. **Divergence and expression excess are dead.** The families with the highest residual `de` in simulated HG00733, where donorsim builds no
   donor-only copies, are LGALS9C, KRT17P1, KRTAP21 and CES1. So this signal measures alleles and paralogy, as register §6.1 found for FARDIV.

## 1. What the frozen tags look like (descriptive, seen, no truth)
| quantity (tied reads unless stated) | HG002 real | HG00733 sim |
|---|---|---|
| tied reads / families | 33,976 / 94 | 12,015 / 174 |
| AS spread 0 (exact AS tie) | 20,710 (61%) | 3,836 (32%) |
| same `subs` on every tied copy (tags carry no copy information; only PSV columns can) | 21,943 (65%) | 4,258 (35%) |
| min-over-copies `subs`, q50 / q90 / q99 | 0 / 3 / 13 | 0 / 4 / 16 |
| min-over-copies `indels` (bases), q50 / q99 | 1 / 20 | 1 / 19 |
| min-over-copies `de`, q50 / q99 | 0.0013 / 0.0101 | 0.0034 / 0.0223 |
| fewest-`subs` copy ≠ best-AS copy (an indel-driven AS preference) | 1,001 | 570 |
| reads with ≥ 1 `foreign` column; recurrent foreign alleles (≥ 2 reads) | 3,002; 554 | 974; 54 |
| **every** tied copy needs ≥ 20 indel bases (a structural residual) | 340 reads / 39 fams (EIF3C 94, SLC25A6 35, ZC3H11A 28, NPIPB2 16) | 119 / 10 (KRTAP21 73) |
| reads with Q ≥ 10 single-copy votes for ≥ 2 copies; clean one-switch (≥ 2 votes each side); recurrent switch intervals | 198; 11; 1 (3 reads) | 102; 6; 2 (5 reads) |
| primary outside the tie (`outside`); MAPQ 0 | 507; 20,873 | 45; 5,062 |
| teacher primaries: clip ≥ 10% of read; with an SA tag; NM per kb q50 / q99 | 1.06%; 1,237; 1.55 / 12.7 | **0.04%; 3**; 3.50 / 20.8 |

**Reading.**
- **Simulation cannot test clip, SA or chimera fingerprints.** donorsim reads almost never clip or split.
- **The structural residual fires in both donors.** NPIPB2's 16 reads match the BRIEF's 57-bp NPIPB2.2 insertion. But EIF3C (r1329's false
  missing-copy family) and KRTAP21 (simulated, with no donor-only copy) fire too.
- **PSV-switch reads are rare.** Both counts are in the table above.

## 2. Which feature separates which class
| class | what its reads look like after O2-style assignment | separating feature | why the obvious tags fail |
|---|---|---|---|
| (i) absorbed from a donor-only copy C* | They tie on, or teach, their nearest copy A. Each read carries C*'s private alleles and structure. | **FP1** (a 3rd haplotype on A); **FP2** when C*'s share is below the caller floor; **FP3 / FP4** when C* differs by an indel or is a hybrid | Above about 10% AF, C*'s alleles enter A's genotype (F2), so the reads stop being foreign or discordant. Mean de or NM mixes C* with A's own alleles. |
| (ii) teachers of a reference copy A the donor lacks | They are a paralog's haplotype and fit A by a few PSVs. Their genotype is learned from themselves, so they show 0 discordance (EXPLORE §5.2). | **FP5**: A is never *anchored*. No read beats every rival by more than the PSVs can explain, and no read crosses A's unique flank or junction. | Every genotype-side test is circular here. G4 is inert and G5 marks the source copy (EXPLORE §5.2). |
| (iii) normal reads | They carry ≤ 2 recurrent haplotypes per copy. Leftover substitutions are near 0 (q50 0). The copy is anchored if its transcript leaves the duplicon. | the absence of FP1-FP5 | — |

## 3. Main fingerprints

### FP1. Haplotype excess on one reference copy (the main fingerprint for (i))
| | |
|---|---|
| (a) signal | Pool copy A's teachers and the tied reads placed on A. Reduce each read to its Q ≥ 10, non-A>G bases at A's VCF records ∪ the family's PSV columns. Group identical vectors into classes of ≥ 2 reads. Two classes conflict when they differ at a shared column. Fire when ≥ 3 pairwise-conflicting classes exist and none equals a tied rival's called genotype on the shared columns. |
| (b) mechanism | A diploid copy carries ≤ 2 haplotypes. A clique of k ≥ 3 conflicting classes proves ≥ k distinct sequences, so some read is not from A's alleles (pigeonhole). Excluding rival-equal classes removes referenced leakage (F1). What remains is a sequence no reference copy explains: an absorbed or collapsed copy (F2). ASE and gene conversion never add a class, so they cause misses, not false fires. |
| (c) literature | Paraphase 2025 §Method step 5: "Total CN = number of unique haplotypes"; per the vault O3 topic note, CN is "not influenced by the CN difference between an individual and the reference". SDA §Method step 3: "repulsion edge for 'mutually exclusive' PSVs". IsoCon §Assumes: "cannot separate alleles … from distinct copies"; k ≥ 3 is the case where alleles stop being an explanation. Tomaszkiewicz: transcript counts "echo the high copy number". Pal 2026: "Sequence isoform diversity … correlated with gene copy number". longcallR §Assumes: diploid, largest phase set kept, so its LowQual unphased 0/1 is only a *symptom* of FP1. |
| (d) confounders | RNA editing (A>G in the transcript frame is excluded). RT template switching, which can make a mixed class recur; add the register's microhomology leg. Leakage from a copy outside the family or the tie (check secondaries). Expressed pseudogene reads. Allele-specific expression (misses only). |
| (e) threshold | Combinatorial: a clique of 3, with ≥ 2 reads per class (the flag's existing constant). It is a certificate, not a score, so its size can rank candidates. |
| (f) compute | Read `teach/<fam>.bam` plus the tied reads' `columns.tsv` rows. The conflict graphs are tiny, so an exact max-clique is enough. A haploid copy (male X/Y non-PAR) needs only k ≥ 2, the same idea as r1337's one-copy-one-group rule. |
| (g) register | Not killed. Adjacent: the ≥ 2-haplotype "raw hidden-copy flag" fired on 7.39% of single-copy genes vs a 7.93% background (§6.1, heterozygosity-dominated); k ≥ 3 is built to exclude that. **r1329 T5:** `--psv-clusters` called a missing copy in 6,907 EIF3C reads (fails at 0.687), so EIF3C must be a pre-stated negative control. G1 (LowQual 0/1) was costly only as an *abstention gate* (EXPLORE §3a). |

### FP2. Teacher discordance excess (DuploMap's rule moved to teachers; (i) at low share)
| | |
|---|---|
| (a) signal | Compute, per teacher of A, its residual alleles: Q ≥ 10 bases outside A's called genotype at A's records and PSV columns. Fire when residual alleles recur (≥ 2 teachers) and are linked (≥ 2 teachers carry both). This is `teachers.tsv: teacher_flag_reads`. |
| (b) mechanism | C*'s reads teach A, and C*'s alleles are visible as discordance only while its share is below longcallR's floor (alt AF ≥ 0.1 and ≥ 2 alt reads). Above that, they become A's genotype (F2 signature). FP2 and FP1 therefore cover complementary share ranges. |
| (c) literature | DuploMap fulltext, Methods: "Reads with high discordance can be result of missing duplicated sequences in the reference genome … Reads that overlap at least five PSVs and show a high rate of discordance are assigned a low mapping quality". Discussion: "does not currently account for missing sequences or copy number changes". |
| (d) confounders | Teacher error ≈ 2% (DM 1.9%, r1330 2.1%) shows up as rival-equal residuals, so exclude those. Editing. Weak hets (SMN2 AF 0.14). Teachers are selected by their alleles (F5). |
| (e) threshold | Combinatorial (recurrence 2, linkage 2). DM's "≥ 5 PSVs" is a hand-set constant; use recurrence instead. |
| (f) compute | Already in `teachers.tsv`. Recompute it with the residual taken against the genotype *and* against every rival's genotype. |
| (g) register | As an O2 copy gate it does not separate bad copies: 0/8 vs 10/79 on HG00733, 7/9 vs 12/21 on HG002 (EXPLORE §2). As an O3 screen it is untested. |

### FP3. Recurrent structural residual, plus rescoring in event units (for (i), and for the outside-the-tie NPIPB2.2 case)
| | |
|---|---|
| (a) signal | Find an indel event from the CIGAR (and MD) that every tied copy needs, at an identical breakpoint and sequence in ≥ 2 reads (EIF3C, NPIPB2 and KRTAP21 in §1 are examples). Then ask whether a secondary alignment just outside the 10-AS window explains the event. A copy wins on events when it has strictly fewer substitutions and no more indel events. If such a copy exists, widen the tie (an O2 repair). If not, flag the event as an unexplained structural allele (O3). |
| (b) mechanism | An insertion in the donor's copy B is charged about 8 AS per indel (memory, corrector) against 5 per mismatch. That pushes B out of the tie, so the read ties on copies that are all worse (NPIPB2.2, 57 bp). A donor-only copy with an SV relative to its nearest paralog leaves the same recurrent event. |
| (c) literature | DuploMap §Methods: a PSV is any differing substring between unique anchors, so indels count. PP25: indel PSVs are polymorphic ("PMS2CL-specific" duplication). longcallR / Clair3-RNA indel precision is 82.7% / 60.5%, so rely on read recurrence, not the caller. Expected failures F8. |
| (d) confounders | HiFi homopolymer error (recurrence plus a non-homopolymer context). Alternative 5′ / 3′ splice sites (NAGNAG): drop events within a few bp of an N operation, where Clair3-RNA's F1 falls from 90.31 to 83.98% at 2-16 bp. Variable TES / poly(A) at read ends: drop terminal events. **The main confounder is a donor SV allele of a present copy:** KRTAP21 73 in a simulation with no absent copy. Only FP1 separates the two, because an SV allele adds no 3rd class. |
| (e) threshold | Event counts with a strict "fewer" are combinatorial. Recurrence is ≥ 2. |
| (f) compute | `reads.tsv` (`indels`, `subs` per copy) gives the screen. The event test needs the input BAM's secondaries outside the window, a light re-read. |
| (g) register | Not killed as O3. `--seq-tie` (r1332) counted unaligned bases, not events. The r1350 substitutions-only arm is an assignment rule and failed on donor substitutions. |

### FP4. PSV-switch (hybrid or converted) reads
| | |
|---|---|
| (a) signal | Order the Q ≥ 10 single-copy votes of `columns.tsv` along the read. Fire on a clean one-switch run (A…A B…B, ≥ 2 votes per side) whose switch interval (between two consecutive PSVs) recurs in ≥ 2 reads. |
| (b) mechanism | A donor-only copy that is a fusion or hybrid of A and B is matched by neither copy alone. A converted tract inside a present copy gives *two* switches instead. |
| (c) literature | Paraphase 2025 §Method step 7: fusions are "haplotypes whose flanks match two different genes; breakpoint = switch at curated PSVs". Parascopy §Assumes: "gene conversion and hybrids break" constant ParalogCN. |
| (d) confounders | RT template switching (register 6.1 M: recurrence alone mis-calls it, so the microhomology leg is required). Het alleles at PSV sites: 42.5% of SD SNVs sit on PSVs with the other copy's allele (DuploMap Key Points). Chimeric reads (SA). |
| (e) threshold | Combinatorial: monotone run structure plus recurrence. |
| (f) compute | `columns.tsv`, done in §1: 11 and 6 one-switch reads, 1 and 2 recurrent intervals. That is low power. |
| (g) register | Fusion category (memory 10-17) and r1347 (SM392 two-copy chimeric primaries) are related. As an O3 flag it is untested. |

### FP5. Anchoring deficit of a copy (the only RNA-only lever on (ii))
| | |
|---|---|
| (a) signal | Copy A is anchored when ≥ 1 read carries evidence only A can produce. Either (1) a read beats every rival by more than the PSVs inside its span can contribute, which means it covers A-unique sequence: a flank, exon, junction or indel PSV. Or (2) a read aligns concordantly across A's duplicon boundary. Fire on copies that are tied and unanchored while their rival is anchored, and in particular on EXPLORE's H5 pairs (separated only by private hom-alt calls). |
| (b) mechanism | A phantom A is taught by a haplotype of paralog D. Its teachers win only by the few PSVs where that haplotype resembles A. Any read that left the duplicon would carry D's flank and go to D. |
| (c) literature | Dishuck §Method step 2: copies validated by "contiguity … plus ≥ 30 kbp unique flanks". Guitart: "structurally equivalent if 90% or more … could be mapped". DuploMap step 2: unique-k-mer LCS, where a single survivor gets MAPQ 254. Greshnova: tied reads "covered transcripts only partially and thus could not capture the small difference". |
| (d) confounders | **Transcripts lying entirely inside the duplicon** (recent duplicates): the copy is *unanchorable*, which is not the same as phantom, so abstain. Low expression (a 0 is weak, so rank it against the rival's anchored fraction). Copy-specific alt-splicing in a real copy (that helps: it anchors). Clipping at the boundary is a symptom, not the test. |
| (e) threshold | Existence of ≥ 1 read is combinatorial. The PSV-only bound is computed exactly from the copy-copy alignment (the per-PSV score swing summed over the span). Nothing is tuned. |
| (f) compute | Teacher AS on A and on every rival, from the input BAM's secondaries (the frozen `reads.tsv` covers tied reads only), plus O1's copy-copy alignment for the boundaries and PSVs. |
| (g) register | Adjacent kill: "the missing-copy signature is CLIPPING" was refuted for (i) (FARCLIP median 0.0006; "no published collapse detector uses clipping"). FP5 is the opposite, a *positive* test for (ii), and is untested. H5 (EXPLORE §5.3) costs present homozygous copies; H5 ∧ ¬anchored is the proposed repair. |

## 4. Considered and not recommended
| feature | (a)/(b) idea | (c)/(d) | (e)/(f) | (g) verdict |
|---|---|---|---|---|
| Per-read residual substitutions (single read) | Edits no candidate genotype explains | An error and a private allele of a present copy look the same | `reads.tsv` subs, `foreign` columns | Only useful when recurrent and linked, which is the existing flag (C1 FAILS after relabel; 107 no_counterpart unflagged) and FP1/FP2. |
| Edit-rate or divergence excess (de, NM) per copy or family | Absorbed reads inflate divergence | Driven by indel error and donor alleles; tops simulated families with no absent copy (§0.3) | rank of median de | **Killed (adjacent):** FARDIV is a paralogy detector, TPR 0.0588 below 0.01 divergence (§6.1). Tied copies are ≥ 99.5% identical. |
| Coverage or expression anomaly (teacher counts, depth ratio) | An absorbed copy raises depth | "expression is not dosage"; 16/104 zero baselines | — | **Killed:** §6.3 TRAP. R2 (least-taught copy) flips sign between donors (EXPLORE §3b). |
| Soft clips at the host | Absorbed reads clip at the flank | Whole-genome alignment gives them a better home | — | **Killed:** §6.1 (FARCLIP 0.0006); 6.2 (no detector uses clipping). Yoo: 0.7% clipped > 200 bp is a mappability statistic. |
| MAPQ or secondary count | Low MAPQ = suspicious | A function of the same AS gap; 61% of HG002 tied reads are exact ties | `calls.tsv` mapq | No new information. DuploMap: 1.9% of MAPQ ≥ 30 reads are wrong anyway. |
| Private or escaping splice junctions | A collapsed copy brings its own junctions | 58.6% of intact nodes carry one | — | **Killed** as a standalone signal (§6.1). Useful only as FP5 anchoring evidence. |

## 5. Evaluation plan (truth used later, only for scoring; nothing tuned here)
- **Which donor tests which class.**
  - HG00733 can test only (ii): FP5, and H5 ∧ ¬FP5. Score them against `build.json`'s "copies without an ortholog" and "silent copies".
  - HG00733 cannot test FP3 (clip/SA) or FP4 (0.04% clips, 3 SA).
  - HG002 tests (i): FP1-FP4 on the families with the 42 donor-insertion reads. It also tests FP3's widening on NPIPB2.2's 16 reads.
  - Both donors are now seen, so their results are descriptive. A fresh HPRC donor gives the verdict.
- **Metrics, all rank-based with no cut.**
  - Precision@k of copies ranked by certificate size (FP1 clique size, FP5 unanchored-vs-rival).
  - Recall of truth copies.
  - Negative controls, pre-stated: EIF3C (r1329), PAR families (allele-rich), KRTAP21 (SV alleles), and the haploid X/Y copies (k ≥ 2 rule).
- **Pre-register as one O3 screen:** FP1 ∪ FP3 for (i) and FP5 for (ii), each with ≥ 2 reads as its only constant. Report FP2 and FP4 alongside, without bars.
