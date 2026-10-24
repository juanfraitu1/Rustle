# Pre-registered expected failures: O2 from the sample's own variant calls (2026-10-22)

**Status.** Committed BEFORE the development run (Task 6) and before any panel run (Task 9), at the user's request (2026-10-22):
"explore [the ThesisVault] and see if we can corroborate what we are doing plus pre-register possible expected failures".

**How it was made.** A blind literature survey by a subagent:
- read the spec, the SDD ledger and memory, all 28 vault paper notes and the cited full texts;
- read no panel output and no truth file.

The controller adopted it unchanged below the line, apart from this header.

**Already addressed when this was committed** (so these predictions are about the code as it will run):
- **§0.7, F2 and F16(3) ("placed `no_counterpart` / `elsewhere` reads are invisible to A1/E1").** Fixed by spec amendment
  c63bbd6a and the Task 5 fix round:
  - every tied read with a definite truth is scored;
  - those placements count as wrong;
  - every arm is scored on the same reads.
- **F16(4): T6's policy choice ignored K1.** T6 will also report K1 for each genotype policy (P1/P2/P3). The choice rule stays
  dev precision, then fixes. A policy whose dev K1 exceeds 2% is excluded.
- **§2.18 reported-only counts.** Adopted as a diagnostics step in Task 6/9, reported only, no new bars.

**Caveat to carry into every claim (§1 row 2, F15).** longcallR's PacBio MAS-seq model was trained on PacBio's HG002 Kinnex
cDNA set. That set is very likely the HG002 reads of panels A and C, so HG002 is not held out for the caller (to verify).

**Novelty framing (§3.3).** "DuploMap's reliable-PSV principle (Prodanov & Bansal 2020, DNA) carried to long-read RNA, with
an assign-or-abstain certificate and an RNA-side missing-copy flag." It is not a new principle. Search for an RNA analogue
before claiming "first".

---

# O2 from the sample's own variant calls: vault survey, pre-registered expected failures, novelty delta

**Date:** 2026-10-22.
**Destination:** to be committed as `docs/PREREG_o2_variant_calls_expected_failures_2026-10-22.md` BEFORE T6 and T9 run.
**Blind:** I read no panel output, no truth file and no directory under `o3_graph/o2vc/`. The SMN2 dev-smoke facts come only from
the ledger and the Task-4 report. These predictions must not be tuned toward.

**What I read:**
- the spec, including "Context: related work and known gaps";
- the SDD ledger, the Task-4 and Task-5 rulings, and the Task-4 report;
- memory `project_o2_variant_calls_2026-10-22` and `reference_paraphase_parascopy`, plus the O2 notes they rest on (r1329, r1330,
  r1333, r1334, r1336, r1343-r1348, the corrector note, `feedback_metric_traps`);
- the docstring of `bench/o3_graph/donorsim.py`;
- every vault paper note (28) and the Topics notes;
- the full text behind every claim cited here.

**Citation key** (short name = vault paper note):
- DM = Prodanov_2020_DuploMap; PSC = Prodanov_2022_Parascopy;
- PP23 = Chen_2023_Paraphase_SMN; PP25 = Chen_2025_Paraphase;
- LR = Huang_2026_LongRNA_SNP_Phasing (longcallR); C3R = Zheng_2025_Clair3-RNA;
- SDA = Vollger_2019_SDA; CH = Chaisson_2017_Multicopy_Duplications_Polyploid_Phasing; IC = Sahlin_2018_IsoCon;
- DV = Shaw_2025_devider; PH = Hosseini_2025_pHapCompass;
- CZ = Canzar_2016_Multimapping_Conflicts; DS = Deschamps_2020_MultiMapped_Review;
- PAL = Pal_2026_Y_Multicopy_Survival; DSH = Dishuck_2025_NPIP; GU = Guitart_2024_TBC1D3;
- GR = Greshnova_2024_YAG_T2T_Isoforms; RH = Rhie_2023_T2T_Y.

`§` names the full-text section. `rNNNN` is a negative-results-register row, quoted from memory. "(inf.)" marks my own inference,
not a paper's claim.

## 0. What the literature adds to the spec (read first)

1. **The core loop is published, for DNA.** DM genotypes reference PSVs per copy from the sample's confidently placed long reads. It
   decides only with PSVs the sample confirms at both copies, and gives MAPQ 0 when no location beats all the others. O2-VC is that
   loop moved to RNA:
   - restricted to AS-tied reads;
   - with a hard certificate and a linked-allele flag;
   - without DM's indel PSVs and without DM's iteration.
2. **Teacher error is about 2%, measured twice.**
   - DM §Results: after PSV re-scoring, 1.9% of HG002 HiFi reads with MAPQ ≥ 30 fall below 10.
   - r1330: 2.1% of non-PAR untied HG002 reads lie off the donor's copy, concentrated in CNV SDs.

   The leaked allele's allele fraction is εr/(1+εr), where ε is the teacher error and r how much more the rival copy is expressed
   (inf., §2.0). With ε = 0.02 it crosses longcallR's 0.10 floor once the rival copy has about 5.6 times the teachers. SMN2's
   AF-0.14 het is the expected product (r ≈ 8).
3. **longcallR is used off-label.** LR filters MAPQ < 10. It trains and benchmarks only inside GIAB high-confidence regions, at
   coverage ≥ 20. Its PacBio MAS-seq model was trained on PacBio's HG002 Kinnex full-length RNA set, very likely panels A and C's own
   reads (verify that `hg002_kinnex` is that PacBio release). C3R's stratified benchmark ranks SDs among the two worst strata, and there longcallR averages about 66% F1.
4. **Caller misses become reference votes.** A missed call or a HomRef record at depth ≥ 6 becomes {reference base}, which is an N3
   vote. In SDs, caller false negatives therefore re-import exactly the reference-PSV errors VC exists to remove. Substitutions alone
   were a coin flip outside SMN/THOC3: 17 fixes, 18 breaks (memory, corrector).
5. **"Dense = collapsed copy" is unsupported for the copies O2 sees.** Reads tie only between copies ≥ 99.5% identical (memory, 10-19).
   A collapsed copy of that kind adds ≤ 0.5 candidates per 100 bp, against the Dense trigger of ≥ 5 (inf.). Dense will mark alignment
   artifacts.
6. **Teachers are selected by their alleles.** A copy's reads that carry its rival's allele are the reads that tie, so they rarely
   teach. DV names the mechanism: reference bias. The genotypes are pushed toward separating alleles (false homozygotes), and a copy's
   second haplotype can look "foreign" (false flags).
7. **Three scoring facts the bars hide.**
   - With fewer than 50 VC placements, A1 ≤ 2% means zero wrong.
   - Placed `no_counterpart` reads are invisible to A1 and E1 (metric trap 09-22).
   - On `outside = 1` reads N1 is wrong by construction, which inflates fixes. T6's criterion (dev precision) ignores K1.

## 1. Corroboration

| # | Element | Supports | Contradicts or qualifies | Does the literature change anything? |
|---|---|---|---|---|
| 1 | Teachers by an AS margin (≥ 10 over every family rival) | GU §Methods "Iso-Seq and transcript analyses", the same rule: reads "with primary minimap2 alignment scores of 10 or greater for a given paralog cluster relative to all other cluster mappings were retained, while other mappings were marked as ambiguous and ignored". DSH §Results: the best mapping must beat the next "by at least one additional mismatch". DM §Methods step 2: PSV genotypes come from reads "aligned with high confidence". | GU and DSH measure the margin against panels of assembled human haplotypes (GU: "all validated human haplotypes"; DSH: 169), not against one reference. A margin against CHM13 copies does not certify the donor's copy. DM §Results: 1.9% of MAPQ-≥ 30 HiFi reads demoted. PP25 §Results: "even reads with high MAPQ can be misaligned due to reference genome artifacts, common CNVs, and high rates of gene conversion". r1330: 2.1% off the donor's copy. r1334: GTF2IRD2's and HSPA6's sample clusters were labelled by near-tie reads, because margin-≥ 10 reads were too few to label them. | Keep the rule, since it is the field's. Treat the teachers as about 2% contaminated, most where expression is skewed and copy number varies (§2 F1-F3, F5). GTF2IRD2's thin margin-≥ 10 support is a known risk for A2. |
| 2 | Per-copy genotypes from an off-the-shelf RNA caller | DM §Methods "Identifying reliable PSVs" makes the same move in DNA: per-copy PSV genotypes from confident reads, and only PSVs that are reference-consistent at both copies are used. LR §Results: PacBio precision 98.5-99.0%, recall 91.0-95.3%, with an A-to-I class. | Off-label on four counts. (i) LR §Methods "candidate SNP sites": MAPQ ≥ 10, secondary and supplementary reads dropped; the design runs `--min-mapq 0`. (ii) LR §Methods "Training dataset": candidates outside GIAB high-confidence regions are removed, and the PacBio MAS model was trained on "the Kinnex full-length RNA dataset (HG002a:MAS:cDNA) from PacBio", so HG002 is not held out for the caller. (iii) LR §Methods "Benchmarking": accuracy is measured only where coverage ≥ 20 inside GIAB high-confidence regions. (iv) C3R §Results "Performance analysis by genomic context": "Low complexity" and "Segmental duplications" have the lowest F1, and longcallR averages about 66% over PacBio's difficult strata (CDS 91.74%, homopolymers 37.07%). | The headline accuracy does not carry over to SD copies at depth 6. Report caller QC on the teachers (§2 F15), and treat HG002 as caller-training data. |
| 3 | Disjoint-allele vote with assign-or-abstain | DM §Methods "Read assignment using PSVs" has the same beats-every-rival structure: the best location b must satisfy P_bi ≥ P_ib for every i, else "we keep the original alignment of the read and assign it a mapping quality of 0". DV §Methods: "If two or more equally good paths exist for a read, we do not assign the read to any path". PH §1: "haplotype mapping uncertainty in the presence of identical haplotype segments". DS §3: the field's defaults are fractional; rescue is "equivalent to a single iteration of the expectation maximization algorithm". A shared allele separating nothing matches DM's constant term for unreliable PSVs and PSC's f ≥ 0.95 rule. | "One separating column is enough" is the most permissive evidence rule in the vault. DM caps each PSV's weight (ŝ = 10⁻³, "to reduce impact of a single PSV on read mapping"). CH (Abstract) assumes "a sizable number of paralog-specific variants". SDA's attraction edges need reads that span ≥ 2 PSVs. samplecopy (r1333) needed ≥ 2 decisive sites. The no-counter-vote veto makes one erased column decisive: on the SMN2 smoke, one AF-0.14 het produced 688 placements. | The shape of the rule is well supported, but its precision rests on single-column genotypes. Report the outcomes of `vc_cols == 1` placements separately (§2 F1, F10). |
| 4 | SNVs only, indels ignored | LR §Results: Clair3-RNA's indel precision is 82.7% (HG004 MAS cDNA) and 60.5% (ONT cDNA), "much lower than the precision of SNP calling". C3R §Discussion: indel accuracy "remains challenging". DM drops low-complexity PSVs, the homopolymer class. | DM §Methods "identifying PSVs": a PSV is any differing substring between unique anchors, so indels are kept unless low-complexity. PP23 §Results: haplotypes are "phased variants (SNPs or indels)". PP25 §Results: a common "PMS2CL-specific" duplication sits on 72.8% of PMS2CL and 31.1% of PMS2 haplotypes, so indel PSVs are polymorphic too. Memory (corrector): AS charges about 8 per indel against 5 per mismatch; HG002's SMN and THOC3 copies differ from CHM13 by a 1-bp insertion; substitutions alone are right on SMN+THOC3 (588/605) but a coin flip elsewhere (17 fixes, 18 breaks). | Yes. Indels create many real ties, so the SNV-only rule both abstains (where copies differ only by indels) and overrides indel-driven AS preferences (§2 F8). DM's high-complexity-indel filter is a ready middle path. |
| 5 | RNA-editing, Dense and low-depth sites treated as unknown | LR §Introduction: editing events "may be mistaken as SNPs". C3R §Results: tagging REDIportal sites raises SNP F1. IC §Discussion: "RNA editing may generate transcripts that have substitutions but originate from the same copy, leading to overestimating the true copy number", which supports excluding A>G from flags. Depth 6 is LR's own candidate floor. Using "unknown" rather than a reference stand-in is more conservative than DM's prior p(0/0) = 0.95. | Dense does not detect collapsed copies at O2's identities. SDA §Results resolves collapses below 0.5% divergence, and O2's tied copies are ≥ 99.5% identical: ≤ 0.5 candidates per 100 bp against Dense's ≥ 5 (inf.). Depth 6 is below every vault benchmark (PP25 §Results, 10× per haplotype; LR, ≥ 20; C3R, zygosity errors fall as depth rises). RnaEdit labels flip between identical runs (ledger, 4/41), and on the smoke 21 of S19SMN2.0's 26 records were RnaEdit. | Keep "unknown". Drop the claim that Dense marks collapsed copies: it will mark artifacts (§2 F2). Expect run-to-run flips (T10). |
| 6 | Weak hets (SMN2, AF 0.14) | Reasons to doubt them: C3R §Methods "Variant zygosity flipping" uses VAF_HET 0.5 with tolerance ε 0.2, so AF 0.14 lies outside the het band (0.3-0.7) and within ε of VAF_Artifact 0.08. LR's candidate floor is alt AF ≥ 0.1 with ≥ 2 alt reads, so 0.14 is barely a candidate. Row 1's ≈ 2% error predicts exactly such AFs (§2.0). | Reasons to believe them: C3R §Results says zygosity flipping is "common in RNA data", through allelic imbalance and isoforms. PP23's "c" haplotypes and DM §Results, where 16.5% of HG001's two-copy PSVs are (0/0, 0/1), make real hets at PSVs common. No vault method separates leakage from allele-specific expression at a single site. | P1 risks placements created by erasure (A1, K1). P3 risks false homozygotes from allele-specific expression (§2 F4). T6's choice should also score K1 (§2 F16). |
| 7 | Foreign-allele missing-copy flag | DM §Methods "DuploMap algorithm": reads discordant at ≥ 5 PSVs get low MAPQ, as possibly from "missing duplicated sequences in the reference genome". In SDA, linking PSVs needs reads that carry both. PP25 Fig. 1: "an extra haplotype … indicates a copy number gain". PSC §Results: OTOA's aggregate copy number is 6 where the reference has 4. SDA: 51.5% of resolved collapses "diverged". | IC §Discussion: "no functionality … that would distinguish allele-specific variants from distinct gene copies". PP23 phased haplotypes into alleles only with pedigrees. SDA tells PSVs from alleles by single-copy depth, which RNA lacks. DV §Results: minority "haplotypes" were mixes "well supported by 'chimeric' reads". DM §Discussion: "does not currently account for missing sequences or copy number changes", so a copy absorbed into a reference copy is never foreign. | The flag means "unexplained haplotype". In a diploid that includes the untaught allele of a present copy (§2 F5), and it misses absorbed copies (§2 F2). C1 is weak evidence in either direction. |
| 8 | No iteration or EM | DS §3: EM and rescue are fractional. Thesis: r579 (EM gained nothing), r1351-1359 (rescue/EM harmful on real reads), r1345/r1346 (echo). | DM §Methods iterates placement and PSV reliability "until the read assignments do not change"; it is DM's only teacher-cleaning step. PSC estimates reliability by EM across samples. | No change. But O2-VC has no way to correct its ≈ 2% wrong teachers. The middle path is one non-fractional teacher check: DM's discordance demotion, applied to teachers. |
| 9 | No facility location | CZ §1: reads are "assigned (uniquely) … such that the total weight of the assignment is maximized", with no abstention. Its 0.19 bound is on the total weight, not on any read's correctness (CZ note §Caveats). Among AS-tied reads the weights tie, so the optimum is indifferent (inf.). Thesis: facility location was 47-59% precise. | Nothing in the vault. | No. |
| 10 | Reference-only null arm (N3) | PSC §Results: when population-unreliable PSVs are present, GATK and Freebayes precision fall from > 0.96 to 0.56/0.59 and recall from 0.52/0.55 to 0.25/0.29. That is the all-reference-PSV against reliable-PSV contrast. DM §Introduction: unreliable PSVs reduce "sensitivity and specificity". r1348: an abstaining metric needs a null arm. | The literature predicts that the sample calls mainly cut N3's breaks (precision). A3 instead tests fixes, the less-supported direction. N3 shares VC's known/unknown gating, so it isolates what the genotypes contain, not the gating. | Report N3's and VC's breaks next to A3 (reported only). |
| 11 | A hand-picked showcase plus a held-out donor | PP25 (259 people, 36 trios, 47 HPRC assemblies), PSC (2,504 + 698 genomes, MLPA, trios) and DM (genome-wide simulation, GIAB) all validate far beyond one locus, so holding out a donor is the right instinct. | Metric trap 09-18: hand-picked sets overstate F by about 0.27. GTF2IRD2 was picked after the corrector (19/0) and samplecopy (31/0) had succeeded on it. Panel E is donorsim output: reads "split evenly over its haplotypes' orthologs", no editing, and no donor copy without a CHM13 ortholog. Its PRAMEF/TRIM43B failure is already known (r1346). DM's real-data gains are checked "only indirectly" (DM note §Caveats). | Panel A shows the method can work, not how accurate it is. Panel E is optimistic for allele-specific expression, editing and missing copies (§2 F13). |

## 2. Pre-registered expected failures

### 2.0 Conventions and derived quantities

**Fields:**
- `calls.tsv`: `vc vc_set vc_cols n3 n3_set flag outside margin mapq`.
- `<family>.columns.tsv`: `kind` ∈ `separating_call`, `separating_ref`, `erased_ref`, `foreign`, `editing`. Per copy it gives copy,
  chrom, pos, ref, genotype and depth, with `.` meaning unknown.
- `<family>.teachers.tsv`: FILTER counts per copy and `teacher_flag_reads`.
- `<family>.flags.tsv`; `vcf/<copy>.vcf` (FILTER, and FORMAT GT, DP, AF, PS).
- `<family>.case.tsv`: `truth_class`, `truth`, `vc_outcome`.
- `<family>.summary.md`: the `vc_cols == 1` count, the evidence class of each fix, and the placements that depend on erased columns.
- `score.<panel>[.<focus>].json`: arms, K1, K2, flags and bars; `run.json`.

**Leakage arithmetic (inf.).**
- Suppose a fraction ε of rival B's confident reads lands on copy A, and B has r times A's teachers at a site. The leaked allele's AF
  on A is then εr/(1+εr).
- It passes LR's floor of 0.10 at εr ≥ 0.111, and the preset's 0.15 (ledger) at εr ≥ 0.176. Two alt reads need ε·n_B ≥ 2.
- With ε = 0.02 (DM 1.9%, r1330 2.1%): r ≥ 5.6 for 0.10, r ≥ 8.8 for 0.15, and n_B ≥ 100.
- SMN2's AF 0.14 at DP 79 (≈ 11 alt reads) fits ε = 0.02 with n_B ≈ 550 and r ≈ 8.

**Small-n A1.** With fewer than 50 placements, wrong/placed ≤ 0.02 means zero wrong. GTF2IRD2 had about 35 decidable tied reads in
the corrector's domain, and minimap2 was wrong on 23 of them.

**Bars:**
- A1: wrong placements ≤ 2% of VC's placements. A2: ≥ 10 fixes. A3: VC fixes more reads than N3.
- K1: ≤ 2% of donor-tied (`own_tie`) reads placed. K2: flagged reads with a present copy ≤ 2%.
- C1: at least one flag holds ≥ 2 `no_counterpart` reads.
- E1: wrong ≤ 2% of placements. E2: fixes > breaks.

### F1. Teacher contamination by a referenced rival (cross-copy leakage)
| | |
|---|---|
| (a) mechanism | About 2% of rival B's confident reads land on A with margin ≥ 10. On a copy expressed 5-9 times less than its rival, they push B's allele past the caller's floor at reference PSVs. The resulting record (usually a LowQual 0/1 whose ALT is B's reference base) makes A = {a, b}, and the column becomes `erased_ref`. Every tied read whose only counter-vote sat there is now placed: on SMN2, one site gave 688 placements. Such a placement is a break when the erased column was the read's correct vote and its remaining votes come from a biased genotype (F5). |
| (b) literature | DM §Results: 1.9% of MAPQ-≥ 30 reads demoted. PP25 §Results: misalignment despite high MAPQ. r1330: 2.1% off the donor's copy, in CNV SDs (GPAT2, SMN2, THOC3, FCGR2A, WASHC1, TBC1D3B). LR §Methods: AF ≥ 0.1 and ≥ 2 alt reads. Arithmetic in §2.0. |
| (c) signature | In `vcf/<copy>.vcf`, a 0/1 record (mostly LowQual) with FORMAT AF 0.10-0.25 on the copy with fewer teachers (`teachers.tsv`), whose ALT equals the paired copy's reference base. In `columns.tsv`, `erased_ref` at that position across most of the family's tied reads. In `summary.md`, "VC placements that depend on erased columns" > 0. In `case.tsv`, breaks over-represented among erased-dependent placements. |
| (d) prediction | Pushes A1 toward FAIL (one break fails GTF2IRD2, §2.0), raises K1, and raises E1 (donorsim's lognormal expression often gives r ≥ 5). Falsifiable: (i) in HG002, at least one non-dev pool family carries such a weak het on its lower-expressed copy; (ii) wherever ≥ 5 breaks exist, the break rate among erased-dependent placements is at least 3 times the rate among other VC placements. |
| (e) likelihood | High: it already happened on SMN2, and ε ≈ 2% is measured twice. |
| (f) mitigation | Partly, through T6's frozen policy: P3 restores the counter-vote, but only below AF 0.15. Follow-up: let a het erase a column only if the copy's own teachers carry the weak allele without the rival's other alleles (a teacher-side linkage test, computed already for `teacher_flag_reads`). Alternatively, demote teachers that disagree with their copy at ≥ 2 separating columns, once, without iterating (DM's discordance step). |

### F2. Teacher contamination by an unreferenced copy (absorption)
| | |
|---|---|
| (a) mechanism | A donor copy C* with no CHM13 counterpart teaches its nearest reference copy A whenever its reads beat every rival there by ≥ 10. C*'s private alleles then enter A's genotype (0/1 at an AF equal to C*'s share, phased as A's "second haplotype" when A is otherwise homozygous). C*'s tied reads therefore carry no foreign allele, are never flagged, and are placed on A or abstain. The scorer cannot see this: a placed `no_counterpart` read is not decidable, so it is in neither the numerator nor the denominator of A1 or E1 (metric trap 09-22). Recent duplicates of a reference copy, the commonest extra copies, are absorbed by construction (inf.). |
| (b) literature | DM §Discussion: "does not currently account for missing sequences or copy number changes". PP25 §Discussion: Paraphase "is not influenced by the CN difference between an individual and the reference", the property O2-VC lacks. PP25 §Results: three CN-variable groups look like false duplications "due to missing genes in the CHM13 T2T assembly". PSC §Results: 15 loci above the reference copy number in more than half the samples (OTOA, 6 vs 4). SDA: 51.5% of resolved collapses diverged. r1343 and r1346: BRANCH/REMOVED worlds and PRAMEF. |
| (c) signature | `case.tsv` rows with `truth_class = no_counterpart` and a non-empty `vc`. In the VCF, ≥ 2 PASS 0/1 records sharing one PS on a copy, at AF 0.15-0.35, with ALTs that are no family copy's reference base (a collapsed haplotype). `teacher_flag_reads` > 0 in `teachers.tsv` only when C*'s share is below the caller's floor. Dense stays at 0-2 records per copy and away from PSVs. |
| (d) prediction | Pushes C1 toward FAIL; A1 and E1 are blind to it by construction. Falsifiable: in HG002, VC places a non-zero number of `no_counterpart` tied reads in at least one pool family. Panel E cannot show F2, because donorsim simulates only orthologs of CHM13 copies. |
| (e) likelihood | High that it occurs; medium for its effect on C1. |
| (f) mitigation | None: it is the spec's "Known limit". Follow-up: a Paraphase-style haplotype count from longcallR's phase sets, where a copy carrying a phased minor haplotype of private alleles becomes unknown and is flagged on the teacher side. Report `no_counterpart` placements. |

### F3. The donor lacks a reference copy, or one donor copy sits under two reference copies (the PRAMEF case)
| | |
|---|---|
| (a) mechanism | The donor's copy D teaches two reference copies A and A′, each read going to the one it resembles more. If D is homogeneous, A and A′ learn the same alleles, the columns become `erased_ref`, and VC abstains, as the spec expects. If D's two haplotypes resemble different reference copies, hap1 teaches A and hap2 teaches A′. One donor copy then yields two disjoint genotypes, a "phantom separation", and D's tied reads are placed by haplotype: half of them land on the copy that is not D's lifted counterpart. |
| (b) literature | r1346: 31 of the corrector's 40 breaks are PRAMEF/TRIM43B; HG00733 lacks PRAMEF8, and the donor's PRAMEF7 reads fill PRAMEF8. r1336 (Y ampliconic genes): two RBMY1 groups share one label, alleles taken for copies. PP23, Table 1 note: the only copy-number error was "two of the three haplotypes being identical in sequence". PP25 §Results: in 22 GRCh38 groups, > 95% of people have copy number 2. PAL: HG002 carries 45 TSPY copies against a median of 32. |
| (c) signature | `separating_call` columns in `columns.tsv` where both copies' genotypes differ from their own reference base, each holding the other's base or a donor allele (complementary swaps). In `case.tsv`, breaks in both directions inside one family. |
| (d) prediction | Pushes E1 toward FAIL, since E contains PRAMEF/TRIM43B; risk to A1 in HG002 in r1330's CNV families. Falsifiable: in E, VC abstains on ≥ 90% of PRAMEF/TRIM43B tied reads and makes ≤ 5 breaks there (corrector v1: 31). |
| (e) likelihood | Medium. |
| (f) mitigation | Erasure covers the homogeneous case; nothing covers teaching split by haplotype. Follow-up: when two copies' genotypes are swapped relative to the reference at ≥ 2 columns, treat them as one copy and abstain. |

### F4. Allele-specific or monoallelic expression
| | |
|---|---|
| (a) mechanism | Take a donor het at copy A, a/b, where b is rival B's base. If b is expressed below the caller's floor, A is called {a}; the column separates A = {a} from B = {b}, and A's minority b-reads vote for B. A read with no other vote for A is placed on B (a break), and own_tie reads from that haplotype get placed (K1). Fully monoallelic a is harmless, since no b-reads exist. Fully monoallelic b makes A = {b} = B, which erases the column (also harmless). P3 (weak alt → {ref}) creates this case on purpose. |
| (b) literature | C3R §Results "Zygosity flipping analysis": "also common in RNA data"; ignoring zygosity raised F1 by 2.43 points (ONT) and 2.88 (PacBio), and removed 4,334 false positives and 4,331 false negatives in BCM Iso-Seq HG004. C3R §Methods: a het "appears as a homozygous reference or a homozygous variant due to allelic imbalance". LR §Introduction: "skewed allele fraction caused by allele-specific expression". |
| (c) signature | A break whose deciding column in `columns.tsv` (often `vc_cols == 1`) has the true copy homozygous, while that copy's teacher BAM carries the read's base in 2-15% of reads. This needs an audit pileup; HomRef records expose FORMAT AF. None in E, because donorsim splits each copy's reads evenly over its haplotypes. |
| (d) prediction | Raises A1 and K1 in HG002; E is unaffected, so part of any precision gap between HG002 and E is allele-specific expression. Falsifiable: if VC breaks at least once in HG002, at least one break is decided at a column where the truth copy's teachers carry the read's base at AF < 0.15. |
| (e) likelihood | Medium: it needs a het at a PSV, strong imbalance and a sole vote. |
| (f) mitigation | None: the GT is taken whatever the AF. Follow-up: no single-site rule separates allele-specific expression from leakage (Row 6). Use linkage among the teachers instead. If the weak allele co-occurs with the rival's alleles at other columns on the same teacher reads, it is leakage; if it is independent of them, it is a real het. |

### F5. Teachers selected by their alleles (reference bias of the teacher set; intrinsic to the design)
| | |
|---|---|
| (a) mechanism | Whether a read teaches depends on its alleles. Copy A's reads that carry B's allele at a PSV lose about two mismatches' worth of margin and tend to tie, so A's teachers under-sample exactly the haplotype (or conversion tract) of A that resembles B. Three effects follow. (i) False homozygotes at donor-het PSVs: F4 without any allele-specific expression. (ii) The private SNVs of A's B-like haplotype are missing from A's genotype, so that haplotype's tied reads carry linked "foreign" alleles and are flagged although A is present (K2). (iii) In a conversion tract, A's teacher depth falls below 6, so A is unknown there and VC abstains. |
| (b) literature | DV §Methods: "rtoa and ator differ owing to reference bias: If a read comes from a haplotype with a true alternate allele, it may systematically align with the reference allele incorrectly". PP25 §Methods realigns every read of a group to one archetype, so no read is chosen by where it was placed. GU and DSH measure margins against panels of assembled haplotypes, so copies are compared against real copies, not one reference. DM re-places reads iteratively. r1329: 73% of HG002 tied reads also tie in the donor. r1334. |
| (c) signature | In `teachers.tsv`, a copy's teacher count far below its share of tied reads. In `flags.tsv`, flags whose reads have `truth_class` unique or own_tie on one copy (the K2 numerator) and whose alleles all sit on one haplotype. Separating columns where about half of a copy's tied reads carry the rival's base: a haplotype, not sequencing error. |
| (d) prediction | Raises K2 in HG002 and E, and lowers A2 through effect (iii). Falsifiable: if any flag holds ≥ 3 reads of a present copy, its alleles will be heterozygous in the donor's copy (checked against the donor assembly on the truth side), not errors. |
| (e) likelihood | Medium-high: the mechanism is structural; its size is unknown. |
| (f) mitigation | None. The disjoint rule confines effect (i) to reads without another vote. Follow-up: genotype each copy from all of its family's reads realigned to one archetype (Paraphase), and use the teachers only to name the haplotypes. |

### F6. RNA editing
| | |
|---|---|
| (a) mechanism | A-to-I editing (A>G on the transcript strand) at an A/G reference PSV makes the A-copy's edited reads look like the G-copy's. It is harmless where longcallR labels the site RnaEdit (unknown) or calls it 0/1 (erased). It is harmful where the editing level in the teachers stays below the floor, so the A-copy = {A}, and an edited tied read has that column as its only vote. Weak T>C and A>G records flip FILTER between identical runs, so a column's known/unknown status, and the placements that depend on it, can flip too. |
| (b) literature | C3R §Results: 8,049 (Iso-Seq) and 15,659 (MAS-Seq) high-quality editing sites; about half of the PacBio sites are intronic; REDIportal tagging raises F1. LR §Introduction. IC §Discussion: editing inflates copy counts. Ledger: 21 of 26 records RnaEdit on S19SMN2.0, and flips in 4 of 41 runs. |
| (c) signature | RnaEdit counts in `teachers.tsv`; `editing` rows in `columns.tsv`; breaks decided by an A↔G column (in the transcript frame) whose A-copy genotype is {A}; differences at RnaEdit/LowQual sites between T10's five re-calls. |
| (d) prediction | At most 1 HG002 break attributable to editing. E has none, since no editing is simulated, so E is optimistic. T10: if S19GTF2IRD2 has any weak A>G or T>C record, at least one VC placement or flag changes across the five re-calls. |
| (e) likelihood | Low for errors; high for power lost at A/G PSVs; medium for run-to-run instability. |
| (f) mitigation | Treating RnaEdit as unknown and excluding editing from flags covers the main cases; editing below the caller's floor is not covered. Follow-up: an A↔G column (transcript frame) may not be a read's only vote. |

### F7. Gene conversion and mosaic copies
| | |
|---|---|
| (a) mechanism | A donor copy converted toward its paralog. A homozygous tract gives it the paralog's base, the column becomes `erased_ref`, and VC abstains or decides elsewhere while N3 votes wrong. An allelic (het) tract is also erased. Mosaic reads vote for both copies and abstain. VC errs here only through F5, when the converted haplotype is the untaught one. |
| (b) literature | PP25 §Results: more than 75% of African PMS2/PMS2CL haplotypes are partly or fully converted, and 36 trios show 4 de novo conversions, 2 of them non-allelic. PP23 §Results: "c" haplotypes. PSC §Discussion: its constant-copy-number assumption breaks under "gene conversion events and hybrid alleles"; VC decides per column and makes no such assumption. DM §Discussion: "gene conversion is well known to result in overlap between PSVs and variants". PAL: 99.94% and 99.87% homogenization. |
| (c) signature | Runs of adjacent `erased_ref` columns on one copy, homozygous for the paralog's base. Fixes whose evidence class in `summary.md` is "reference PSV corrected by the sample". `n3` wrong on the same reads where `vc` is right or empty. |
| (d) prediction | Pushes A3 toward PASS. In E, where HG00733's conversions are real, N3 makes at least twice as many breaks as VC. |
| (e) likelihood | High that conversion occurs; low that it causes VC errors. |
| (f) mitigation | Mitigated by design (`erased_ref`); what remains goes through F5. |

### F8. Real indel PSVs lost (SNV-only)
| | |
|---|---|
| (a) mechanism | Many real-read ties are made by an indel: AS charges about 8 per 1-bp indel and 5 per mismatch. VC sees only SNVs. Where the copies differ only by indels inside the common window, it abstains. Elsewhere it overrides an AS preference that came from an indel, on SNV evidence. That is right when the indel is a donor-vs-reference difference (SMN, THOC3) and wrong when it is a real PSV the donor carries. N4, which uses left-aligned indels, can decide reads VC cannot. |
| (b) literature | LR §Results: SNVs only; Clair3-RNA's indel precision is 82.7% and 60.5%. C3R §Results "Indel performance": PacBio indel F1 averages 87.46%. DM §Methods: PSVs of any length, with low-complexity ones filtered out. PP23: "SNPs or indels". PP25: the PMS2CL duplication on 72.8% and 31.1% of haplotypes. Memory (corrector): 588/605 right on SMN+THOC3, 17 fixes and 18 breaks elsewhere. |
| (c) signature | `reads.tsv` rows whose tied copies differ in `indels` but not in `substitutions`, giving an empty `vc` and `vc_cols` 0. N4 fixes, in `score.json`, on reads VC abstains on. Among VC placements that disagree with `minimap2`, the placed copy's alignment carries more indel bases, which follows from how the tie arises. VC breaks concentrated in families whose CHM13 copies differ by a non-homopolymer indel inside the common window. |
| (d) prediction | Lowers A2. On GTF2IRD2, VC fixes fewer reads than N4's 19; falsifiable: VC makes ≤ 15 fixes there. The sign of E2 is unaffected. A1 and E1 rise only in families with indel PSVs. |
| (e) likelihood | High for lost power; medium for errors. |
| (f) mitigation | Not mitigated: it is known gap 1. Follow-up: DM-style high-complexity indel PSVs as columns under the same disjoint rule, with their own prereg. |

### F9. Low expression of a copy
| | |
|---|---|
| (a) mechanism | A tied copy whose teacher depth is below 6 at a column is unknown there, and an unknown copy can never be beaten. Every read whose tie set contains it abstains. The foreign-allele test needs every tied copy known, so the same copy also blocks flags. Many-copy families (NBPF1, NPIP, SPDYE1) are blocked by their least-expressed tied copy, which is also where leakage AF peaks (F1). |
| (b) literature | LR §Methods: a depth floor of 6, and a benchmark at ≥ 20. PP25 §Results: 10× per haplotype. C3R §Results: DP ≥ 8 instead of ≥ 4 adds 4.19 F1 points (ONT). IC §Results: DAZ reaches > 90% recall only at about 410×. DSH §Results: 55.6% of reads individually assignable, even with 169 haplotypes. r1334. |
| (c) signature | An empty `vc` in `calls.tsv` whose `vc_set` contains a copy with genotype `.` at every column of the read in `columns.tsv`. That copy's count in `teachers.tsv`; `failed_copies` in `run.json`. |
| (d) prediction | Pushes A2 toward FAIL: on GTF2IRD2 I guess 3-9 fixes, and on NBPF1 A2 needs 10 of about 11 minimap2-wrong reads. Pushes C1 toward FAIL: SPDYE1 is blocked, and NPIPB2 needs both .5 and .7 known. Some families in E place no read at all. |
| (e) likelihood | High. |
| (f) mitigation | Correctly conservative; there is no fix inside the rule. Follow-up: before scoring, report per family the tied reads whose tie set is fully known at ≥ 1 column, so an A2 miss is classed as a power failure, not an error. |

### F10. 5′-truncated or short reads
| | |
|---|---|
| (a) mechanism | A column must be aligned in every tied alignment, and separating sites cluster in a few exons (for SMN, exons 7-8 and downstream). Truncated or short reads, or reads whose tied alignments clip or splice differently, cover no separating column and abstain. A read that covers a single column becomes a `vc_cols == 1` placement, resting on one genotype and one base. |
| (b) literature | GR §Discussion: multi-mapping reads "covered transcripts only partially and thus could not capture the small difference between gene copies". IC §Discussion: "impossible to separate transcripts originating from copies with identical exonic sequences". PP23 §Results, for exons 1-6: "read alignments are often ambiguous in this region, even for long reads". C3R: F1 drops from 90.31% to 83.98% within 2-16 bp of splice boundaries (Iso-Seq). |
| (c) signature | Abstentions dominated by `vc_cols = 0`; VC's placement rate rises with `ms` (`reads.tsv`); the `vc_cols == 1` share in `summary.md`. |
| (d) prediction | Lowers A2; precision is unaffected apart from single-column errors. Falsifiable: in pooled HG002, VC's placement rate among tied reads in the top `ms` tertile is at least twice the rate in the bottom tertile. |
| (e) likelihood | High for lost power; low for errors. |
| (f) mitigation | Correct abstention. The exon-only footprint is RNA's built-in ceiling compared with DNA methods. |

### F11. Near-identical copies (DAZ)
| | |
|---|---|
| (a) mechanism | DAZ copies are near-identical, and DAZ's 2.4-kb repeat units lie inside one gene, so tied alignments can pair a read base with non-homologous repeat positions (spurious columns). Most DAZ reads tie on all copies, so teachers are scarce, copies are unknown, and VC abstains. The Y is haploid: a 0/1 or 1/2 record on a non-PAR Y copy cannot be a germline het, so every het on DAZ marks a copy mixture, a collapse or a repeat misalignment. A119b's Y is not CHM13's (CHM13 v2.0's chrY is HG002's T2T-Y), so A119b's AZFc copies can differ in number and structure. |
| (b) literature | PAL §Results: palindrome arms 99.94% identical, arrays 99.87%; DAZ has 2-19 repeats. RH: the T2T-Y is "from the HG002 genome", and the DAZ "entire repeat array, consisting of 2.4 kb composite units … falls within one gene annotation". GR §Discussion: "the DAZ gene includes almost identical repeats that could cause alignment issues". PP25: Y groups indistinguishable by sequence alone. IC: DAZ needs about 410× for > 90% recall. r1336: 393 of 579 DAZ tied reads tie on all four copies, and DAZ2 and DAZ3 differ at one expressed site. |
| (c) signature | Panel B's `score` (`--mode none`): placed, abstained, flagged, and `vc_cols == 1`. Deciding positions in `columns.tsv`; 0/1 and 1/2 records in the DAZ copies' VCFs; `flags.tsv`. |
| (d) prediction | B is descriptive. ≥ 95% of DAZ tied reads abstain. Every VC placement has `vc_cols == 1`, and the deciding columns fall at ≤ 2 sites (expected: r1336's DAZ2/DAZ3 site). N3 places at least as many reads as VC. At least one DAZ copy has a het record. A DAZ flag whose alleles lie inside the repeat units is an artifact. |
| (e) likelihood | High for the abstention; medium for a het record. |
| (f) mitigation | Mitigated by abstaining. Haploidy is a free control for contamination: report het records on the non-PAR X and Y copies of male samples. |

### F12. Reference errors or minority alleles in CHM13 at PSVs
| | |
|---|---|
| (a) mechanism | A CHM13 base error at a PSV, or CHM13 carrying a minority allele: the copy's teachers carry the population base, the copy is called 1/1, and the column is erased (VC right, N3 wrong). A CHM13 false duplication: reads tie exactly, no teacher exists, the copy is unknown, and VC abstains. A copy missing from CHM13 is F2. |
| (b) literature | DM §Methods: "some PSVs may be assembly artifacts". PP25 §Results: 22 GRCh38 false duplications, most already flagged by CHM13, and three groups where CHM13 lacks a gene. PP23 §Results: GRCh38 carries the minority SMN1 haplotype S1-2. Yoo_2025 §Key Points: read-depth conclusions about copy number overturned by complete assemblies. |
| (c) signature | 1/1 records whose ALT equals the paralog's reference base at the same CHM13 site in HG002, HG00733 and A119b wherever they are covered: a non-reference homozygote shared across independent donors. |
| (d) prediction | Helps A3 and E1 and harms no bar; such sites will exist in at least one pool family. |
| (e) likelihood | Low that it causes harm. |
| (f) mitigation | Mitigated (`erased_ref`). |

### F13. Simulation artifacts in HG00733 (panel E)
| | |
|---|---|
| (a) mechanism | donorsim splits each copy's reads evenly over its haplotypes (no allele-specific expression). It simulates no editing and no library chimeras. It drops copies whose orthologs overlap on a haplotype, and exons that change by > 50 bp. It never simulates a donor copy that lacks a CHM13 ortholog. So E cannot show F2, F4 or F6, and shows only part of F3. Its truth is exact, so a placed read from donor-identical copies is right or wrong by chance. E is also not pristine: r1346 ran the corrector on it, and the PRAMEF/TRIM43B failure shaped the spec's "Known limit". |
| (b) literature | The `donorsim.py` docstring, read as the method's definition. PP25 §Results: the simulation floor (10 kb, 10×, 0.05%) shows how simulations set favourable conditions. DM note §Caveats: real-data gains checked "only indirectly". CH note §Caveats: simulation only, with the sizable-PSV assumption true by construction. DV note §Caveats: accuracy simulated. |
| (c) signature | E's 0/1 records sit at AF ≈ 0.4-0.6, and RnaEdit records appear only at simulated errors. E has fewer erased-dependent placements per placement than HG002. Every flag in E counts toward K2. |
| (d) prediction | E's wrong fraction is below HG002's. Falsifiable: if both have ≥ 100 VC placements, HG002's wrong fraction exceeds E's. E says nothing about C1. |
| (e) likelihood | High. |
| (f) mitigation | Not mitigated. Follow-up: re-simulate with beta-binomial allelic ratios, editing at REDIportal sites, and donor-extra copies. |

### F14. Truth-label artifacts (liftover)
| | |
|---|---|
| (a) mechanism | HG002's truth lifts the donor's copies to CHM13 through a chain. In SDs and CNVs a chain can fail: a gap gives `no_counterpart` although a CHM13 copy exists, so C1 counts a spurious hit and K2 is understated. A chain can also pair a donor copy with the paralogous CHM13 copy: the truth is swapped, so a right VC placement scores as a break and a wrong primary as right. chrM has no chain and is excluded. CHM13 v2.0's non-PAR Y is HG002's own, so any `no_counterpart` there is a lift artifact. |
| (b) literature | DSH note §Method in brief: paralogs are named by phylogeny, "not … synteny alone", because interlocus gene conversion breaks positional identity. GU: minigraph could not group TBC1D3 copies. PP25 §Results: most disagreements with HPRC assemblies were traced to assembly errors. RH: the T2T-Y comes from HG002. r1330: 258k chrM `no_counterpart` reads. r1348: the FCGR2A lift was never verified. |
| (c) signature | Breaks whose read aligns with fewer edits to the donor locus of VC's copy than to the truth copy's (an audit re-alignment, on the truth side). C1 flags whose donor locus aligns to a CHM13 copy (C2). Breaks in both directions inside one family. |
| (d) prediction | Spurious breaks in A1, spurious passes in C1, and an understated K2. Falsifiable: at least one of any HG002 panel-A breaks turns out to be a lift artifact on audit. |
| (e) likelihood | Medium. |
| (f) mitigation | Partly mitigated by C2 and the chrM exclusion. Follow-up: before reporting, re-align every break, and the reads of every counted flag, to both donor haplotypes. |

### F15. Caller false negatives and off-label behaviour
| | |
|---|---|
| (a) mechanism | Where a copy's teachers carry a non-reference allele but longcallR writes HomRef, or the allele misses the 0.10-AF / 2-read floor, the copy is "known = {reference base}". That is an N3 vote, so caller misses in SDs re-import the reference-PSV errors VC exists to remove. Other quirks: non-determinism (4/41); unphased hets written as LowQual, whose GT is still used; approximate phasing; `--min-mapq 0`, which the authors never benchmarked; and HG002 Kinnex in the training set. |
| (b) literature | LR §Results: recall 91.0-95.3%, and only in high-confidence regions; it "may miss true SNPs due to phasing difficulties". LR §Methods: the phasing "may miss the global minimum"; trained on HG002a:MAS:cDNA, in high-confidence regions only. C3R §Results: SDs among the lowest-F1 strata, longcallR averaging about 66%. Memory (corrector): substitutions alone go 17/18 outside SMN/THOC3. |
| (c) signature | HomRef records with FORMAT AF ≥ 0.10 on a copy; `separating_ref` columns at those sites; breaks decided there; differences between T10's re-calls. |
| (d) prediction | Raises A1, E1 and K1. Falsifiable: at least one HG002 pool family has a HomRef record with AF ≥ 0.3 at a reference PSV. |
| (e) likelihood | Medium. |
| (f) mitigation | Not mitigated. Follow-up: treat a HomRef record whose AF reaches the caller's own candidate floor as unknown; this uses the caller's constant, not a new one. |

### F16. Scoring artifacts
| | |
|---|---|
| (a) mechanism | (1) With fewer than 50 placements, A1 means zero wrong. (2) On `outside = 1` reads N1 is wrong by construction, so every right VC placement counts as a fix; NPIPB2's 10 reads show such reads exist. (3) Placed `no_counterpart`, `poor_match` and `elsewhere` reads fall outside A1 and E1 (a universe intersection). (4) T6 picks the genotype policy by dev precision on decidable reads, breaking ties by fixes. own_tie placements (K1) are not in its criterion, and P1 maximizes the placements created by erased counter-votes (688 on SMN2). (5) If SMN2 and MSTRG3150 count among "HG002's pool families", data used for T6 enters K1, K2 and C. (6) K1's denominator, about 9.6k non-PAR own_tie reads (r1329), means about 190 placements fail it, and one SMN2-like family can supply them. |
| (b) literature | Metric traps 09-18 (hand-picked sets) and 09-22 (a universe intersection inflates the worst arm's precision); r1348 (null arm). |
| (c) signature | `arms.VC.placed` in `score.A.S19GTF2IRD2.json`; fixes split by `outside`; the truth classes of placed reads in `case.tsv`; K1 per policy in the T6 doc. |
| (d) prediction | A1 is a knife-edge, `outside = 1` inflates A2, and one family can decide K1. |
| (e) likelihood | Certain: these are arithmetic. |
| (f) mitigation | The reported-only counts in §2.18. |

### 2.17 Panel verdicts: best guesses, to be checked and never tuned toward

P is my probability that every bar of the panel passes.

- **A, GTF2IRD2: FAIL**, most likely on A2. I guess 3-9 fixes: the teachers are too thin at depth ≥ 6 (r1334), and VC sees only SNVs where
  N4 made 19 fixes. A1 holds at 0 wrong, and A3 passes if VC places anything. P ≈ 0.15 with K1 and K2 included; 0.3 for A1-A3 alone.
- **A backup, PDPR: the best of the three**, P ≈ 0.35. In r1334, margin-taught labels gave 37/0, and minimap2 is wrong on 34 reads. A1
  is a zero-error knife-edge.
- **A backup, NBPF1: FAIL, underpowered.** A2 needs 10 of about 11 minimap2-wrong reads, and its many copies leave unknown copies that
  block placements.
- **K1 (HG002 pooled): PASS**, P ≈ 0.7. It fails only through erased-counter-vote placements in an SMN2-like family under P1 or P2.
- **K2 (HG002 pooled): PASS**, P ≈ 0.85. The SMN2 smoke flagged 2 of 1,655 tied reads, and the denominator is large.
- **B, DAZ: as expected.**
  - ≥ 95% of tied reads abstain.
  - At most a few `vc_cols == 1` placements, at the one DAZ2/DAZ3 site.
  - At least one het record on a DAZ copy.
- **C-human, C1: FAIL**, P(pass) ≈ 0.4. Absorption, the every-copy-known rule and the linkage rule each remove candidates. A pass would
  come from GTF2H2C or NPIPB2, and must clear C2.
- **E, HG00733: lean PASS**, P ≈ 0.5. E2 ≈ 0.9 and K2 ≈ 0.9; E1 ≈ 0.6, failing (if at all) in families with skewed expression and in
  PRAMEF/TRIM43B.
- **D and C-gorilla: not run** (the fast path), so this round has no one-genome no-harm check.

### 2.18 Reported-only counts that make these predictions checkable (no new bars)
1. VC placements of `no_counterpart`, `poor_match` and `elsewhere` reads (F2, F16).
2. Per family, before scoring: the tied reads whose tie set is fully known at ≥ 1 column (F9).
3. VC's and N3's breaks next to A3 (Row 10).
4. Fixes split by `outside` (F16).
5. Per copy: 0/1 records with AF < 0.25, and whether each ALT is the paired copy's reference base (F1); HomRef records with AF ≥ 0.10
   (F15).
6. Het records on haploid copies in male samples (F11).
7. K1 on the dev families under each T6 policy (F16).
8. Right and wrong among erased-dependent placements, wherever a truth exists (F1).

## 3. Missing prior art and the novelty delta

### 3.1 Side by side

| | DuploMap (DM) | Paraphase (PP23, PP25) | Parascopy (PSC) | O2-VC |
|---|---|---|---|---|
| Substrate | Long DNA reads (CLR, HiFi, ONT) in SDs > 5 kb at ≥ 97% | HiFi DNA, 160 GRCh38 groups > 10 kb at > 99% | Short-read WGS, 167 loci, cohorts | Long RNA reads (Kinnex FLNC), AS-tied reads of named CHM13 families |
| What it decides | Location and MAPQ of every read whose primary overlaps an SD | Haplotypes (= copies), copy number, phased variants, fusions | Aggregate and per-paralog copy number per sample; each PSV's f | Each AS-tied read's copy, an abstention set, or a flag |
| Evidence | Reference PSVs of any length (high-complexity), pair-HMM windows, unique-k-mer LCS | Every variant site of all reads realigned to one archetype | Read depth plus PSV allele counts | SNV columns common to all tied alignments |
| Whose variants | The sample's: per-copy PSV genotypes from reads with MAPQ ≥ t; reliable = posterior (0/0, 0/0) | The sample's own haplotypes | The population's: f across ≥ 50 samples (f ≥ 0.95) | The sample's: per-copy GT from margin-≥ 10 teachers via longcallR |
| Sample-private alleles | Not used (reference PSVs only) | They define the haplotypes | Not used | Can separate (hom-alt) or erase (het) |
| Decision rule | Likelihood; the best location beats each rival pairwise; per-PSV cap | Graph assembly, then read → haplotype (tie handling not described) | HMM, then EM, with quality scores | Combinatorial: assign c iff c beats every rival (≥ 1 vote, no counter-vote) |
| Abstention | MAPQ 0, original alignment kept | Groups of low diversity left inseparable; no-calls | Low quality or no per-paralog estimate | An explicit set of unbeaten copies |
| Certificate | None (a probability) | Haplotype plus a diagnostic base | None | Deciding columns listed per placement |
| Missing copy | Reads discordant at ≥ 5 PSVs demoted; copy-number changes "not accounted for" | An extra haplotype = a gain; immune to the reference's copy number | Aggregate copy number above the reference (OTOA) | ≥ 2 linked recurrent foreign alleles in tied reads |
| Iteration | Yes, until assignments stop changing | Not described | EM | None |
| Indels | Yes (high-complexity) | Yes | Not described | No |
| Validation | Simulations (precision 0.995), GIAB HG002 calls | 259 people, 36 trios, 47 HPRC assemblies | 3,202 1kGP genomes, MLPA, trios | One hand-picked HG002 family; one simulated donor (HG00733) |

### 3.2 What overlaps, what O2-VC adds and what it lacks

**DuploMap (the closest analogue).**
- Overlap:
  - it post-processes minimap2;
  - it genotypes per-copy PSVs from the sample's confidently placed reads (DM's teachers);
  - only PSVs the sample confirms at both copies decide;
  - a location must beat every rival pairwise, else MAPQ 0;
  - reads discordant at ≥ 5 PSVs are demoted as possibly coming from a duplicated sequence the reference lacks.
- O2-VC adds:
  - RNA: spliced, exon-only, expression-skewed and edited reads;
  - a scope of AS-tied reads only;
  - a stock RNA caller in place of a built-in PSV genotyper;
  - sample-private alleles that can separate, and shared alleles that erase (DM only rates reference PSVs as reliable or not);
  - a hard certificate that lists the deciding columns;
  - a linked foreign-allele flag (≥ 2 recurrent alleles, ≥ 2 reads) in place of a per-read discordance count;
  - a reference-only null (N3) and a selective baseline (N2) on the same reads.
- O2-VC lacks:
  - indel and multi-base PSVs;
  - DM's iteration, which is its teacher-cleaning loop;
  - a per-PSV influence cap;
  - the unique-k-mer LCS pre-filter for reads that reach past the SD;
  - a probabilistic MAPQ;
  - a genome-wide evaluation;
  - any copy-number handling (DM lacks it too).

**Paraphase.**
- Overlap:
  - the evidence is the sample's own copies, not reference positions;
  - copies are treated as haplotypes;
  - it abstains where copies are as similar as alleles;
  - gene conversion is visible per haplotype.
- O2-VC adds:
  - RNA;
  - per-read placement onto a named reference copy. The names come from where minimap2 put the teachers. PP25 §Results: "It is often
    not easy to assign haplotypes of a paralog group to individual genes without prior knowledge", so Paraphase names haplotypes by
    unique flanks, curated PSVs or c.840;
  - a per-read certificate;
  - the tied-read scope.
- O2-VC lacks:
  - phasing: multi-site haplotypes from all reads realigned to one archetype, immune both to placement errors and to copy-number
    differences;
  - copy number from the haplotype count plus depth (in RNA, depth is expression);
  - detection of identical haplotypes, which needs depth (PP23's only copy-number error);
  - fusion calls;
  - population and trio validation.

**Parascopy.**
- Overlap:
  - a reference PSV is a hypothesis;
  - only reliable PSVs count;
  - no estimate is given where PSVs are unreliable (C4A/B: 7 of 50 reliable).
- O2-VC adds: one sample at a time, per-copy genotypes in place of population frequencies, read assignment, and RNA.
- O2-VC lacks:
  - population-level reliability, so it cannot tell a donor het from leakage (Row 6);
  - copy number from depth, which finds missing copies such as OTOA directly;
  - replicate and Mendelian validation.

### 3.3 Novelty, stated honestly

**Not new.**
- Genotyping PSVs per copy from the sample's confident long reads, and deciding only with the PSVs the sample confirms (DM: DNA,
  2020).
- Abstaining when copies cannot be told apart (DM, PP25, DV).

**New, as far as the vault shows.**
- Doing it for long RNA reads.
- With an off-the-shelf RNA caller.
- With a threshold-free combinatorial certificate restricted to AS-tied reads.
- Letting sample-private alleles separate copies.
- With a linked foreign-allele flag.
- With an evaluation design that holds a reference-only null, the selective baseline and own_tie calibration.

The claim to make is "DuploMap's reliable-PSV principle carried to long-read RNA, with a certificate and an RNA-side O3 flag". It is
not a new principle. The vault holds no RNA analogue of DuploMap, but a targeted search must come before any "first" claim.

### 3.4 Methods from the vault to adopt, compare against or cite

**Adopt:**
- DM's high-complexity indel PSVs (gap 1).
- DM's discordance demotion, applied to teachers once (F1, F2).
- A Paraphase-style haplotype count from longcallR's phase sets, to mark collapsed copies (gap 2).
- PSC-style "no estimate" reporting: each family's power ceiling, before scoring (F9).

**Compare against (truth side or reported only):**
- Paraphase run on HG002 HiFi DNA, if GTF2IRD2, PDPR or NBPF1 are among its 160 GRCh38 groups, to cross-check copy counts. It must
  never become an input.
- Clair3-RNA as a second caller, if it installs without Docker; C3R reports it ahead of longcallR in SDs.
- DSH's best-plus-one-mismatch rule needs the donor assembly as input (an anti-goal), so it can only be cited.

**Cite:**
- GU, the source of the 10-AS rule.
- PAL §Limitations: transcriptome and assembly from different individuals prevent "a direct one-to-one transcript-to-copy mapping";
  O2-VC's same-sample premise answers this.
- GR §Discussion (own work): multi-mapping reads left unresolved.
- Tomaszkiewicz_2023_YAG_Isoforms: asks for such tools.
- IC: alleles versus copies, and editing.
- DV: abstaining on ties, and reference bias.
- CZ and DS as foils.
- LR and C3R for the caller's limits.
- Tang_2020_FLAIR: FLAIR's reassignment keeps MAPQ ≥ 1, the same selective baseline as N2.

### 3.5 What the vault lacks
- **Allele-specific mapping bias.** WASP-style correction, and RNA alignment to personal genomes or pangenomes. This is the theory
  behind F5, and the alternative the anti-goals rule out.
- **How common allele-specific and monoallelic expression is** in lymphoblastoid lines and testis, needed to size F4.
- **RNA editing** in SDs and in testis (REDIportal itself, ADAR).
- **Kinnex/MAS-seq library artifacts**: concatemer splitting, chimeras, RT template switching. These are sources of foreign alleles.
- **The minimap2 paper.** The definitions of AS, MAPQ and secondary hits (`-N`, `-p`) that "tied" rests on.
- **How accurate liftover and chains are in SDs**, for the truth side.
- **Gene-conversion rates** (Dumont 2015; Vollger 2023). The vault cites them only through DM and PP25.
- **DM's and Paraphase's supplements**: the discordance filter, the haplotype graph, and how ties are handled.
