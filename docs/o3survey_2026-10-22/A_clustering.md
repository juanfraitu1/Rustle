# A. Clustering reads by PSV pattern: what a missing or extra copy leaves in RNA reads (O3 survey, 2026-10-22)

**Scope.** This file covers the vault notes and full texts of SDA, Chaisson 2017, IsoCon, devider, pHapCompass, Paraphase 2023/2025 and
DuploMap. It also covers the brief's docs (O2-VC Results; `O2O3_EXPLORE` §5; expected failures F1-F5, F7, Row 7), the register rows
listed in §4, and the memory notes on IsoCon, samplecopy and psv-clusters.
**One post hoc look at the frozen outputs**, with no truth read: `columns.tsv` of the tied reads. The scripts are in the scratchpad,
`o3A/hap2.py`.

## 0. Bottom line
1. **Clustering gives only a lower bound on copy number, and that bound can be proved.**
   - Take k reads that pairwise conflict at linked columns. They need k haplotypes, so the donor needs ⌈k/2⌉ diploid copies (k on a male's
     X/Y outside the PARs).
   - Without DNA depth, no clustering result gives an upper bound.
   - So **case (i), an extra donor copy, can be certified** whenever it adds a haplotype beyond what its absorbing copy can carry.
   - **Case (ii), a CHM13 copy the donor lacks, can never be certified by haplotype counts.** It can be certified only by a read that
     links the "phantom" haplotype to the surviving copy (§2, Lemma 2).
2. **Each method in the vault decides copy number in a way that does not survive the move to RNA** (§1).
   - SDA and Paraphase use single-copy DNA depth.
   - Paraphase 2023 and Chaisson's DMP need the ploidy given.
   - Under (ii), every method in the vault calls the phantom a present copy, or would no-call it.
3. **The frozen outputs, read post hoc on tied reads at separating columns only.**
   - **Edges from single-column conflicts give false haplotypes.** On PRAMEF7, a 3-read minor allele produces a 4-clique. HG00733's
     PRAMEF7 reads come from one donor copy (`O2O3_EXPLORE` §0.2).
   - **Edges that need ≥ 2 conflicting columns:**
     - the family-level clique bound exceeds 2N in 0/131 HG00733 and 0/50 HG002 families;
     - PRAMEF7 has clique 2 (one diploid copy suffices);
     - HG002's real insertion families are not exceeded: SPDYE1 2 against N = 4, RASA4B 2 against 2.
     - The family-level excess test is therefore specific but blind. The per-copy test needs the teacher reads' alleles, which the frozen
       outputs do not hold.

## 1. Each method on the two scenarios
(i) = a donor copy C* that CHM13 lacks, whose reads are absorbed by its nearest reference copy A. (ii) = a CHM13 copy A′ that the donor
lacks, taught by the donor's copy D of paralog A (PRAMEF7.1: D's haplotype h2 resembles .1).

| method (vault) | clustering object; how copy number is decided | (i) absorbed donor-only copy | (ii) phantom A′ (PRAMEF7.1) |
|---|---|---|---|
| **SDA** (Vollger 2019 §Methods "PSV graph construction") | PSV nodes; attraction = reads carrying both PSVs; repulsion when "<10% of the mean coverage" carry both; correlation clustering "does not require the number of clusters". A PSV is accepted only if its depth is "consistent with the read depth of unique regions". | DNA: excess depth marks the collapse at A, and the PSVs then split C* off as a "diverged" contig (51.5% of resolved contigs). RNA has no unique-region depth, so PSVs cannot be told from alleles (expected failures, Row 7). | No collapse, so SDA is silent. Clustering A ∪ A′ reads gives 2 clusters, which SDA reads as 2 paralogs: the phantom is accepted. SDA's substrates CHM1/CHM13 are haploid, where clusters are copies; in a diploid, alleles inflate the count. |
| **Chaisson 2017** §2-3 | §2: "a read conflict as two overlapping reads that are non-gap and disagreeing at a site … grouping all conflict-free reads". CC discovers the cluster number. DMP needs S, since MEC "decreases monotonically with ploidy". | CC finds 3-4 clusters on A when the copies are PSV-rich. It assumes "a sizable number of paralog-specific variants"; a recent C* has few. | 2 clusters = 2 copies, so the phantom is accepted. DMP with S = 2 fits both readings (one het copy, or two hom copies). |
| **IsoCon** (Sahlin 2018 §Discussion) | Reference-free nearest-neighbour clusters plus a statistical test; copies = cliques of transcripts that differ only by large indels. | C* becomes its own transcript if it differs by ≥ 1 base at enough depth. Our chain flags it (r1196: 15/16 NPIPA2) and links "alleles" by δ (r1199). Real reads: it detects but does not recover (r1363). | h1 and h2 become 2 transcripts. The δ-link calls h2 ".1's allele", so .1 is reported present. "no functionality … that would distinguish allele-specific variants from distinct gene copies". |
| **devider** (Shaw 2025 §Methods, Results) | Positional de Bruijn graph over (site, allele); it "naturally splits if enough variation is present and collapses under ambiguity". Tied paths go to no haplotype. | An extra haplotype on A if C* differs at the given SNPs. Minority haplotypes can be "a mix of the two majority haplotypes … supported by 'chimeric' reads". The O3 arm was REFUTED (commit 6994b228, r1354): it is reference-bound, and a second copy folds into a mixed haplotype. | Per copy: 1 haplotype each; on the union: 2. Consistent with both readings, and devider has no notion of copy number. |
| **pHapCompass** (Hosseini 2025 §3) | K and genotypes are given. A SNP pair admits only O(K) valid phasings (3 of 256 for a tetraploid). | With K = 2, ≥ 3 recurrent two-site patterns on A contradict K, but the method cannot raise K. | The union reads fit K = 2 exactly (one copy), and K = 4 with identical pairs fits too: "haplotype mapping uncertainty in the presence of identical haplotype segments" (PH §1, as quoted in expected failures, Row 3). |
| **Paraphase 2023** (Chen 2023, Table note) | Haplotypes from all SMN1/2 reads realigned to SMN1, named by c.840. Its one CN error: "two of the three haplotypes being identical in sequence". Its no-call: depth "could not reliably distinguish CN1 vs. CN2 when only one haplotype was found". | An extra haplotype raises CN, unless C* is identical to an allele of A, which collapses it. | Two haplotypes plus DNA depth give CN 1, so A′ is absent. Without depth (RNA), this is exactly Paraphase's own no-call condition. Haplotypes were phased into alleles only with pedigrees. |
| **Paraphase 2025** (Fig. 1; §Discussion) | "an extra haplotype … indicates a copy number gain in the paralog". Total CN = "number of unique haplotypes, adjusted by depth". "not influenced by the CN difference between an individual and the reference". | The design case. In RNA it counts only expressed, distinct haplotypes. | The haplotype count (2) does not show the loss: "adjusted by depth" does. Genes are told apart only where haplotypes reach non-homologous flanks (5 groups; 23 low-diversity groups are indistinguishable). |
| **DuploMap** (Prodanov 2020 §Methods) | No clusters. Each copy's PSVs are genotyped from confident reads; "reliable" = (0/0, 0/0). Reads discordant at ≥ 5 PSVs get low MAPQ, as possibly "missing duplicated sequences in the reference genome". | Seen only if C* is discordant at ≥ 5 reliable PSVs. A recent duplicate agrees with A and stays at high MAPQ. "does not currently account for … copy number changes". | Both copies' PSVs look reliable (each copy is taught by its own haplotype), so h2 reads are placed confidently on A′. This is O2-VC's error, made in DNA (O2-VC is this loop moved to RNA: expected failures, §0.1). |

## 2. The combinatorial core (what can be proved; my synthesis, with Chaisson §2 as the base definition)
- **Conflict graph.** G(R) has one vertex per read. Two reads are joined when they disagree at shared columns. Q ≥ 10 bases only; A>G in the
  transcript frame is never counted, as in VC.
- **Lemma 1** (exact; error-free; one allele per read and column). R can be explained by ≤ m haplotypes iff χ(G) ≤ m.
  - Proof: reads that agree pairwise agree on every column, so each colour class is one haplotype.
  - **The m = 2 case is decidable in polynomial time.** Run a bipartiteness test; an odd cycle is the certificate.
  - **Any k-clique certifies ≥ k haplotypes.**
- **Corollaries for (i).** Each is a certificate of a donor copy number above the reference's.
  - On one copy's reads: an odd cycle in a diploid, or any edge on a male's X/Y outside the PARs.
  - At family level: χ > 2N.
- **Lemma 2: (ii) is non-identifiable.** Take a transcript lying wholly inside the duplicated segment. Compare two donors:
  - Donor 1: copy D carries h1/h2, expressed at allelic ratio α, and A′ is absent.
  - Donor 2: A ~ h1/h1 and A′ ~ h2/h2, with expressions in the ratio α : (1-α).
  - Both produce the same distribution of read sequences, hence the same CHM13 alignments, so any read-only rule gives both the same
    output.
  - **The mirror case is also non-identifiable.** "A heterozygous" and "A homozygous plus a homozygous C*" (case i) are equivalent in the
    same way. So the commonest absorbed copy, a recent homozygous duplicate, is invisible to counting too (F2: "phased as A's second
    haplotype").
  - **What breaks the symmetry, and none of it is in the reads:**
    - DNA depth, which is Paraphase's "adjusted by depth";
    - a molecule that reaches copy-specific sequence (Paraphase 2025: unique flanks; SDA: "22% of SDA contigs transition into unique
      sequence");
    - a prior on α (no allele-specific expression) or on how often copies are homozygous. These are priors, not proofs.
- **Error handling.** Each allele must be recurrent (≥ 2 reads), and an edge needs ≥ 2 conflicting columns.
  - These are the field's minimum evidence units: SDA's attraction edges need both PSVs on a read, VC's flag needs linkage, and samplecopy
    needed ≥ 2 sites.
  - The frozen PRAMEF7 shows why the second condition is needed (§0.3).

## 3. Fingerprints computable after O2-style assignment
| | **FP1 haplotype excess** (case i) | **FP2 unexplained divergent cluster** (case i) | **FP3 pair parsimony** (case ii, screen) | **FP4 anchored ownership** (case ii, certificate) |
|---|---|---|---|---|
| (a) signal | Among one copy's reads (teachers ∪ VC placements), an odd cycle or k ≥ 3 clique at ≥ 2-column conflicts, where the extra class matches no haplotype taught to another copy. Haploid X/Y: any edge. Family: χ > 2N. | A read cluster on a copy whose consensus is farther than allelic divergence δ from every reference copy (`copy_assign --psv-clusters`; IsoCon chain). | Over the reads of a tied pair (A, A′), G is bipartite, each copy's teachers lie in one class, and the two classes differ. One diploid copy explains both: CN_min = 1 < 2. This is the read-level form of H5. | A read carrying class h's linked alleles plus a copy-specific anchor (junction, exon, UTR or end found only at A, or reads past the SD into A's flank) assigns h to A. Flag when every class teaching A′ is anchored to A; clear A′ when any class is anchored to A′. |
| (b) mechanism | (i): C*'s reads teach A, or are placed on it, so A carries ≥ 3 haplotypes. | (i) with C* diverged beyond allele level (SDA's "diverged" class). | (ii): D's two haplotypes teach A and A′ separately (PRAMEF7 per `O2O3_EXPLORE` §5.2). | (ii) when the transcript reaches copy-specific sequence. |
| (c) literature | Paraphase 2025 Fig. 1; Chaisson §2; SDA's CC; pHapCompass's O(K) valid phasings. | SDA "diverged" 51.5%; DuploMap's ≥ 5 discordant PSVs; IsoCon. | Paraphase 2023's no-call (one haplotype, CN1 vs CN2 by depth only); IsoCon §Discussion; Lemma 2. | Paraphase 2025 (genes told apart only via non-homologous flanks; haplotypes sharing reads are grouped into alleles); SDA (22% reach unique sequence). |
| (d) confounders | F1 leakage brings B's haplotype in (excluded by the match rule); chimeras and template switching give mosaics (devider), so the extra class must carry an allele neither other class has; editing; gene conversion makes a real C* look like a mosaic and be missed. | Hyperdivergent alleles (register §6 "Copy-vs-allele can be settled on RNA alone"; §6 MHC row); immunoglobulin (§6 IG rows); δ is an allele prior; reads that are 5′-spread fail to cluster (r1363). | **Present homozygous pairs** (identical by Lemma 2). Allele-specific expression. F5 reduces false flags here, since tied reads add A's third class; this is a falsifiable gain over H5. Converted copies. | Readthroughs and chimeras (r852); isoforms annotated only on a paralog (r1299 CROSS); conversion of the anchor; liftover only matters for truth. |
| (e) threshold-free? | Combinatorial. Constants: ≥ 2 reads and ≥ 2 columns. | No: δ (p99 allele divergence, r1199). | Combinatorial, with the same 2/2 constants. Recall-oriented; a deficit is never a proof. | Combinatorial; one anchored read is a certificate (≥ 2 for robustness). |
| (f) from our outputs | `columns.tsv` covers tied reads only. **The teachers' alleles must be added:** re-run longcallR without `--no-bam-output` (HP tags) or pile up each copy's teacher reads; then BFS for bipartiteness. | Already in Rust (`--psv-clusters`, opt-in). Run it on HG002's pool BAM. | Same inputs as FP1, per tied pair; the H5 pairs (13 HG00733 / 2 HG002) are the first set to check. | Per-read junction/exon sets from the BAM, against the copy-specific features of the O1 catalog. `reads.tsv` lacks them. |
| (g) register | Not killed. r1194 (male chrX false-2 2/439) supports the haploid form; r1196 says an AF floor loses the swamped native haplotype, so count reads, not calls. | Simulation PASS (r1292/r1293). On real reads: detection yes, recovery no (r1363); realignment PASS (r1364). Register §6 "Copy-vs-allele can be settled on RNA alone" (killed: "needs DNA") is the boundary. | H5 is hindsight (§5.3). G4 is inert, G5 harmful, R2 flips sign (§3b). | Untried. Nearest rows: r731 (canonical splice sites do not discriminate) and r1299. |

**Not proposed: allelic balance** (the two classes in a 1:1 ratio means alleles). It needs a binomial α. donorsim's even split flatters it
(F13). Clair3-RNA reports zygosity flipping as "common in RNA data" (expected failures, F4). At most it could be a ranking feature, never a
gate.

## 4. Clustering already tried (register and memory)
| attempt (row) | result | why it failed or held | lesson for clustering |
|---|---|---|---|
| RNA allele count (r1194) | Pilot PASS: male chrX false-2 2/439; single-copy precision 588/590, recall .76. Multi-copy run never done. | "RNA cannot separate homozygous from hemizygous". | The haploid form of FP1 is calibrated. Lemma 2 is this boundary. |
| IsoCon on the family net (r1195; r1363) | Simulation: NPIP 44/44. Real: 13/43 (r1195). r1363 FAIL: 66/94 target reads unclustered; ±1 A poly(A) consensus error. | IsoCon does not converge when 5′ ends spread (C spans 44 kb); its 0.999 class is fragile. | Cluster on PSV columns, not whole-read consensus. Detection (one candidate per unit) did work. |
| Deletion test (r1196) | No fake allele: 97% of the absorbed copy's reads swamp the native copy (native base 2-6%, under the .20 floor). | Expression asymmetry. | (i) can appear as the native haplotype vanishing, not as a third one. Count classes of ≥ 2 reads, without an AF floor. |
| IsoCon linking by δ (r1199, r1200) | Held-out HELP, but any copy closer than δ to a paralog is called "an allele by construction". | δ is an allele prior. | FP2 inherits δ. FP1 and FP3 avoid it. |
| samplecopy EM (r1333; r1334; r1336; r1337) | FAIL 82 right / 4 wrong, at weak donor margins (1-3 sites); purity does not explain the errors. YAG calls are reproducible but are NOT provenance: two groups share a label, and one group covers two copies. | Clusters are labelled by majority, and a tied pair is never separated. | On a haploid chromosome, one copy = one group is a visible failure test. That is FP1's haploid form. |
| devider arm (commit 6994b228; r1354) | REFUTED on 6/6 clauses. | Reference-bound: in 30/53 deletion families, 0% of the missing copy's reads aligned. One contig per haplotype. Mixed haplotypes. | Use a PDBG only for phasing within a cluster. |
| O2-VC flag, G4/G5/H5 (`O2O3_EXPLORE` §5) | C1 is a liftover artifact. G4: 0 pairs. G5 marks the source copy. H5 catches PRAMEF7.1 but costs 66 / 108 right placements. | The flag is blind to absorbed copies (F2) and fires on a present copy's untaught allele (F5). | FP3 is H5 rebuilt on reads. Pre-register it before any look. |
| psv-clusters (r1292, r1293) | Simulation PASS (425/426; 0/152 controls); exact ties alone found 80/310. | Never run on real reads. | Cluster every read O2 places, not just ties. |

## 5. Is case (ii) detectable from RNA reads? (honest judgment)
- **Not in general.**
  - When every expressed molecule lies inside the duplicated segment, Lemma 2 makes a missing copy and two homozygous present copies
    observationally equivalent.
  - The register reached the same boundary empirically: register §6 "Copy-vs-allele can be settled on RNA alone" ("needs DNA"), and r1194 (homozygous vs hemizygous).
- **Signals that exist:**
  - FP4, a molecule that ties a phantom-teaching haplotype to the surviving copy. This is the only RNA certificate.
  - FP3, a parsimony screen: "these reads need one copy fewer than the reference". Its false positives are exactly the present homozygous
    pairs.
- **For PRAMEF7.1:**
  - FP3 fires. The tied-read clique is 2, H5 marks the pair, and .0 has no 0/1 record (explore §5.3).
  - FP4 is probably unavailable, since .0 carries no h2 teacher and so no anchored h2 read exists. That is checkable from the BAM.
  - The honest output is therefore **"possible missing copy; RNA cannot decide"**, which is the ranked screening flag Canzar asked for.

## 6. Evaluation plan (truth used only for scoring; pre-register first)
- **FP3 / FP4 recall: HG00733 `donorsim/.../build.json`.**
  - The positives are pairs where one copy is "without an ortholog" or "silent" and its tied rival has an ortholog.
  - False positives are flagged pairs where both copies have orthologs.
  - Report how many positive pairs have any anchored read: that is FP4's ceiling.
- **FP1 / FP2 recall: HG002's ~42 donor-insertion reads** (SPDYE1 28, RASA4B 8, USP17L11 3, POTEM 2, BMS1P1 1).
  - Exclude NPIPB2's chain artifact.
  - Specificity: HG00733, which has no donor-only copies, so every FP1 flag there is false; and HG002's male chrX outside the PARs (haploid
    FP1).
- **Bars.** Only "≥ 2 reads, ≥ 2 columns", fixed now.
  - HG00733 is spent for H5, so FP3 needs a fresh donorsim donor with beta-binomial allelic ratios (explore §5.4).
