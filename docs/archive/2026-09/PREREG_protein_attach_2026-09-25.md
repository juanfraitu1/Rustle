# Pre-registration: protein attachment, the manual extra-sensitive step of the de novo families

**Written 2026-09-25, before any number of this step exists.** Tool: `tools/protein_attach.py`. Supplementary figure:
`figures/fig_protein_supp.py` (`captions/figS_protein.md`).

User decision (2026-09-25 16:00): *"Protein is NOT in the default. It is an OPTIONAL, MANUALLY INVOKED extra step
('extra-sensitive mode'): for families with disconnected / strongly diverged members, align them at the protein
level to attach missing members. Figures: protein only in a supplementary figure; main figures use external
references."* The default de novo family definition stays the driver's `families` stage (reads -> seeded assembly
loci -> one representative per locus, its read-derived exon coordinates = the positional exon sum -> `mcl_families
--from-gtf`, exon-sum >= 0.60, MCL 2.8). This step reads that stage's outputs and never modifies them.

## 0. What I have already seen (exposure)

- The development inputs' structure, not any outcome: human chr16 de novo families (the Fig. 7 development run,
  `${work}/fig7/current/human_chr16.denovo.fam.*`): 2,802 loci, 111 families with 353 member spans (after
  fold-within-clusters), every locus strand `+` or `-`.
- The genome FASTAs are soft-masked (a 1-Mb window: CHM13 chr16 62% lower case, mGorGor1 36%, mPanTro3 50%).
- Earlier protein arms and their numbers: r1027 (de novo ORFs, development chr16 +2.6 points pair coverage,
  cross-family edges 33 -> 38; held-out annotated-node arms +0.6/+0.0/+0.4), r1028/r1029 (protein-only edges are
  90-97.5% unverifiable, domain/fold sharers), r1030 (coverage floors 0.30-0.90 do not change it), r1031 (the edges
  are real, the truth was missing), r1096 (protein tier precision 0.29-0.67 vs Compara).
- No number of this step, on any substrate. Human chr6 (the held-out verdict contig below) has never been used for a
  protein decision; it was scored once for Fig. 7 (de novo vs guided, no decision).

## 1. How this differs from the refuted protein arms (register check)

| row | refuted | this step |
|---|---|---|
| r1027 | protein edges from translated de novo locus representatives, UNIONED into the family graph before clustering | never a graph edge: no family is re-clustered, no two RNA families are merged; only loci the RNA rule left out can be added, one at a time, to one existing family |
| r1028/r1029 | §6ko protein edges (e <= 1e-5, >= 0.30 of the longer) as a sensitive mode; 90-97.5% of added edges have no nucleotide alignment, mostly fold sharers | the same base hit, but an attachment also needs identity at least as high as the family's own loosest member (section 2.5), a scope floor at the RNA divergence cliff, and a unique passing family |
| r1030 | raising the COVERAGE floor to 0.70-0.90 | coverage stays at §6ko's 0.30 (r1030: no floor separates fold sharers); the discriminator tested here is IDENTITY calibrated per family, which no row has tested |
| r1031 | "the protein edges are false merges" | it asked for a gene-tree truth: this step is judged against Ensembl Compara (gene trees) and Liftoff copies, never against Soto (SD-scoped) or protein-homology families (circular), which appear only as a labelled secondary reference |
| r1027 confound | held-out arms had to translate ANNOTATED nodes (no de novo assemblies off chr16) | de novo families now exist off chr16 (Fig. 7 runs on chr6, chr2, chr8, chr10, gorilla chr10/chr20; genome-wide per sample once the `families` stage has run) |
| T19 / r831 | best-hit counting | best hit only CHOOSES the one family an attachment may go to; every evaluation counts qualifying pairs |
| r637 | label a node by the gene with the most shared bases (large-gene attractor) | a locus takes a gene label only when that gene's exons cover >= 0.50 of the locus's exonic bases; otherwise it is unlabelled (uninformative, never a disagreement) |

## 2. The step (exact definition)

### 2.1 Inputs

The `families` stage outputs of one run, `P = PREFIX.fam`: `P.clusters.tsv` (families, member spans),
`P.loci.gff3` (every de novo locus: its span and its representative's exons), `P.loci.tsv` (loci folded into a
member span), `P.loci.paf` (the all-vs-all; optional, only for the nucleotide-evidence label), and the genome FASTA.

A locus is a **member** of family F when its span is a member span of F or `P.loci.tsv` folds it into one.
Every other locus is **unattached**: loci with no alignment to anything, and loci whose alignments to members failed
the RNA rule (sub-threshold, "weakly connected") or were split off by MCL.

### 2.2 The protein of a locus (no annotation CDS in the de novo mode)

Spliced sequence = the genome bases at the representative's exons (the positional exon sum), reverse-complemented
on the minus strand. Protein = the **longest stop-to-stop open reading frame on the transcribed strand** (three
frames; six only for a locus without a strand), no start codon required, translated with the standard code, an
`N` codon ends a frame, **among the frames whose bases are less than half soft-masked**, kept if it has
**>= 100 amino acids**.

Justification:
- *Transcribed strand only.* IsoSeq reads are stranded full-length cDNA and every assembled locus carries a strand
  (chr16: all 2,802). A protein is encoded only on the transcribed strand; the antisense frames only add chance
  ORFs.
- *Stop-to-stop, no ATG.* This is the definition of the repo's existing de novo protein tier
  (`denovo_pipeline.rs::longest_orf_aa`). The step targets diverged members and pseudogene copies, which lose start
  codons; and de novo 5' ends are dispersed (§6w3: |d5| 136-229 bp median) while the CDS is complete (§6p4: 99.8%),
  so an ATG rule would cut the CDS start of any 5'-short representative.
- *Soft-masked majority excluded.* A frame that is mostly transposable-element sequence (for example an L1 ORF in a
  3' UTR) gives TE homology, the protein analogue of r1098's repeat-glued loci. Majority (0.50) is the only
  parameter-free cut. Masking protocols differ between genomes, so the excluded count is reported per species.
- *>= 100 aa.* The floor of the only earlier de novo test (r1027). Its codon-shuffled null gave 0 edges on four
  chromosomes at this floor, so chance ORFs are not a false-edge source there; a lower floor has no null evidence.
  Length does not certify coding (r1027: every locus has a >= 50 aa frame), so chance is controlled by the E-value
  and the null arm (2.7), not by the length. Loci whose best frame is 50-99 aa are counted and reported as lost.

### 2.3 Search

BLASTP (BLAST+ 2.17.0, the binary `bench/truth.py` uses), exactly `truth.py`'s command: `-evalue 1e-5
-max_target_seqs 100000`, the same 9-column tabular output. **Database = the proteins of every member locus of every
family; queries = the proteins of every locus** (members give the calibration and the merge report, unattached loci
the candidates). Queries run in contiguous shards, each written through a `.tmp` + rename, keyed by a manifest
(md5 of both FASTAs, BLASTP version, arguments): resumable, `--budget-s` stops cleanly and exits 75.

### 2.4 Base hit (the §6ko edge style)

For an ordered pair (query locus q, database locus s), q != s, whose representatives share no genomic base (their
exons do not overlap; a nested gene in an intron is allowed): HSPs chosen greedily by
bit score, non-overlapping on the LONGER protein (`truth.edges_from`), coverage = union of the chosen intervals /
length of the longer protein, identity = sum of identities / sum of alignment lengths over the chosen HSPs,
score = sum of their bit scores. A **base hit** needs coverage >= 0.30 (§6ko's constant, unchanged; r1030 found no
floor better). E <= 1e-5 is applied by BLASTP.

### 2.5 Family calibration (no new constant)

For each member m of family F with a protein: its nearest-sibling identity s(m) = the highest identity of a base hit
between m and another member of F (either direction). F is **calibrated** when at least two members have s(m); its
calibration identity **I_F = min s(m)**: how far F's loosest member already sits from its nearest sibling. A
family whose members are not protein-homologous to one another is not calibrated and cannot receive attachments
(a family defined by shared UTR or non-coding sequence cannot be extended by protein homology).

### 2.6 Attachment rule (assign or abstain)

For an unattached locus u with a protein whose representative exons overlap no member's representative exons (else
verdict `overlaps_member`: it shares genomic bases with a member, so a protein hit would be trivial):

1. F* = the family of u's best BLASTP hit (highest summed bit score over all hits to members). No hit: `no_hit`.
2. F* not calibrated: `family_uncalibrated`.
3. No base hit (coverage >= 0.30) to a member of F*: `below_coverage`.
4. Best base-hit identity to F* below I_F*: `below_family_identity`.
5. That identity below **0.60**: `below_scope_floor`. 0.60 is the measured divergence cliff of RNA edges (§6t1:
   below protein identity 0.60 RNAs do not align) and the band where protein-only homology is dominated by fold
   sharers (r1096: protein-tier precision 0.29-0.67); it bounds the step to the scope of the thesis.
6. Another family also passes steps 2-5: `ambiguous` (not attached; reported as a merge hint).
7. Otherwise `attached` to F*, with its best member, identity, coverage, score, and a label: **nucleotide-corroborated**
   when `P.loci.paf` has any record between u and a member of F* (a sub-threshold alignment), else **protein-only**.

### 2.7 Null arm (every run)

Each candidate query protein is shuffled (a fixed per-locus seed; the amino-acid permutation that a codon shuffle
of the ORF induces) and searched against the same database; the same rule, with the real calibrations, is applied.
Null attachments estimate chance attachments.

### 2.8 Proposed merges (reported, never applied)

Two member loci of different families A and B whose base hit reaches identity >= max(0.60, I_A, I_B) (each I only
when that family is calibrated), and every `ambiguous` locus, are written to `.merges.tsv` with the family pair,
the supporting pairs and the best identity. No family is ever merged (r1028: protein homology links domain and fold
sharers; merging two RNA families is out of this step's remit).

### 2.9 Outputs (separate files; the RNA families are never touched)

`OUT.orfs.tsv` (every locus: family, strand, exonic length, best frame, aa, masked fraction, status),
`OUT.calibration.tsv`, `OUT.candidates.tsv` (every unattached locus with a protein, its verdict and evidence),
`OUT.attached.tsv`, `OUT.clusters.tsv` (the RNA clusters' rows unchanged, plus one row per attached locus, extra
column `added_by` = rna / protein), `OUT.merges.tsv`, `OUT.null.tsv`, `OUT.params.tsv`, `OUT.proteins.faa`,
`OUT.members.faa`, the BLASTP shard directories. (Added while implementing, before any number: `OUT.candidate_hits.tsv`,
every candidate x family it hits, which section 3.5's ceiling needs.)

## 3. Evaluation (supplementary figure only)

### 3.1 Gene labels

A locus takes the RefSeq gene or pseudogene record (the species' `genes.tsv` annotation cache of Figs 6-7) whose
merged exons cover the most of the locus representative's exonic bases, if that is >= 0.50 of them; otherwise the
locus is unlabelled.

### 3.2 References

- **Primary (human): Ensembl Compara release 116, paralogues at any duplication age** (connected components of
  the Compara pairs at every level up to Opisthokonta, over protein-coding CHM13 genes; the `compara_families`
  construction of Fig. 7 at the deepest level). A gene tree separates an ancient paralogue from a fold sharer,
  the distinction r1031 said was missing.
- **Secondary (human): Compara at Primates** (Fig. 7's headline level; the thesis scope of recent families).
- **All species: Liftoff copies** (the Fig. 8 self-lift, `-copies -sc 0.95`): two loci are Liftoff-related when each
  is >= 0.50 covered (exonic bases) by a placement of the same source gene. Liftoff only lifts copies >= 95%
  identical, so it CERTIFIES attachments (a lower bound on precision) and never falsifies one.
- **Secondary, labelled circular: protein-homology families** (Figs 6-7's annotation protein families). They are
  built from protein homology, so agreement with them is partly by construction (r1029, r1031); never a verdict.

### 3.3 Classes of an attachment u -> F (per reference)

- `same_gene`: u's gene is the gene of a member of F (a second locus of a gene F already holds; excluded from
  precision).
- judgeable: u's gene and the gene of at least one member of F are in the reference's universe.
- `true`: judgeable and the reference relates u's gene to a member's gene; `false`: judgeable and it relates none.
- Precision = true / (true + false) (judgeable attachments), with a Wilson 95% interval.

### 3.4 Baseline: the RNA families' own members, same test

Member-level precision of the RNA families on the same substrate and reference: over members whose gene is in the
universe and whose family holds another member with a different gene in the universe, the fraction related to at
least one of them. The attachments must be at least as precise as the members the RNA rule itself admitted.

### 3.5 Sensitivity to missing members

Missing members = unattached loci (overlapping no member, gene in the universe, gene not held by any family) whose
gene the reference relates to the gene of a member of some family F. Sensitivity = those attached to such an F /
all of them. Ceilings reported: missing members with a protein, and with a base hit to that F.

## 4. Substrates and exposure

- **Development (rules may be inspected, not tuned): human chr16.** Every constant above is fixed now; chr16 is only
  checked for crashes, the null and sanity.
- **Held-out verdict (human, now): chr6**, de novo families of the Fig. 7 run (the whole genome-wide A119b assembly
  restricted to chr6). Never used for any protein decision.
- **Descriptive (no verdict until Liftoff exists): gorilla chr10 (NC_073234.2, untouched) and chr20 (NC_073244.2,
  development of the seeding rule).**
- **Genome-wide (when the driver's `families` stage has run per sample):** the held-out samples of the closed loop,
  human_testis (Compara, verdict repeated with chr16 dropped), gorilla_KB3781, chimp_PTR, orangutan_PPY (Liftoff
  certified fraction, descriptive), plus human_A119b and gorilla_OR6737. Species are never pooled.

## 5. Bars (held-out, human, Compara at any duplication age)

| outcome | verdict |
|---|---|
| null attachments > 0 (or > 1% of attachments when there are >= 100) on any substrate | ⛔ **VOID** (chance attachments) |
| Wilson 95% upper bound of attachment precision < baseline member precision | ⛔ **UNSAFE**: the step dilutes the families |
| precision < baseline, interval contains the baseline | ⚠ **INCONCLUSIVE** (too few judgeable attachments) |
| precision >= baseline, >= 5 true attachments in >= 2 families, sensitivity to missing members >= 0.10 | ⭐ **USEFUL**: keep as the manual extra-sensitive step, recommend it for families with diverged members |
| precision >= baseline, fewer true attachments / families or sensitivity < 0.10 | ⚠ **SAFE, SMALL**: keep as a manual option, document the ceiling |
| 0 judgeable attachments | ⚠ **NO EVIDENCE** |

Reported without a bar: the corroborated / protein-only split, the proposed merges and their Compara status,
Primates-level precision, Liftoff certified fractions, the protein-homology secondary.

**Predicted before looking: ⚠ SAFE, SMALL.** The family calibration asks for identity close to the members' own, so
few loci pass; those that pass should be real paralogues (precision at or above baseline). Sensitivity to missing
members should stay below 0.10: r1026/r1101 put most reference pairs outside any alignment of expressed loci, and
r1027's de novo gain was +2.6 points on development.

I will not change the ORF rule, the base hit, the calibration, the 0.60 floor, the references, the classes, the
substrates or the bars after seeing any number. A different rule needs its own pre-registration.

## Outcome (2026-09-25): ⚠ **SAFE, SMALL**, as predicted, on ONE judgeable attachment

Runs: `tools/protein_attach.py` on the Fig. 7 de novo families (`${work}/fig7/current/<species>_<contig>.denovo.fam`),
outputs in `/mnt/linuxdisk/tmp/rustle_figures_dev/protein_attach/`; scoring `figures/fig_protein_supp.py data`
(tables `figures/data/figS_protein_{substrates,refs}.tsv`). No constant was changed after any number. Each run took
1-26 s and < 0.2 GB.

| substrate | loci | families (calibrated) | loci with a protein | candidates | hit a member | attached | null |
|---|---|---|---|---|---|---|---|
| human chr16 (development) | 2,802 | 111 (51) | 1,746 | 1,470 | 28 | **0** | 0 |
| **human chr6 (held-out verdict)** | 5,012 | 54 (11) | 2,456 | 2,339 | 84 | **1** | **0** |
| gorilla chr20 NC_073244.2 (descriptive) | 1,142 | 29 (24) | 1,087 | 1,007 | 213 | 1 | 0 |
| gorilla chr10 NC_073234.2 (descriptive) | 956 | 6 (1) | 895 | 887 | 1 | 0 | 0 |

**Development (chr16).** No crash; null 0; 0 attachments. The 28 candidates that hit a member stop at: below the
family's loosest-member identity 9, coverage < 0.30 9, family uncalibrated 10. chr16's families are near-identical
copies (I_F 0.99-1.00 in most calibrated families), so a candidate at 0.99 identity is refused (DN_chr16_80122322_16:
0.9895 vs I_F 0.9938). Two member-pair merge proposals (Compara: 0 judgeable).

**Verdict (chr6, Compara at any duplication age).**
- Null: 0 attachments on every substrate -> not VOID.
- Baseline member precision 14 of 14 = 1.00 (95% CI 0.78-1.00). Attachment precision 1 of 1 = 1.00 (0.21-1.00): the
  HLA-E locus attached to the HLA class I family (HLA-A/B/C/F/G/K; identity 0.766 over 89% of the longer protein,
  I_F 0.750; a sub-threshold nucleotide record exists, so it is "nucleotide-corroborated").
- True attachments 1 (< 5) in 1 family (< 2); missing members 30, attached 1: sensitivity 0.033 (< 0.10).
- => **⚠ SAFE, SMALL.** Stated plainly: n = 1 judgeable attachment; this is a scope result, not a precision estimate.
- Secondary: Compara within primates, members 8 of 14 (0.57), the attachment 0 of 1 (HLA-E / class I is older than
  primates); no missing member at that level. Protein-homology families (circular): members 14 of 18, attachment 1 of 1.

**Where the step stops (descriptive, observed after the verdict).** Of chr6's 30 missing members, 29 have a protein
and 27 have a counting hit (coverage >= 0.30) to the family Compara relates them to, which is in every case the family
of their best hit. The **family calibration** stops them: 20 are less identical than the family's loosest member
(HLA class II DRA/DQA1/DQB1/DMB/DPA1/DPB1 to the DRB1/DRB3 family, I_F 0.893; MICA/MICB/ULBP1-3/HFE to class I,
I_F 0.750; BTN2A1/2A2/3A1; GSTA3/4; TUBB/TUBE1; CD109) and 6 hit families whose members share no protein homology
(SLC22A2/3/16, ZNF311/391, POPDC1). Most of these are paralogues older than primates, i.e. outside the families' own
divergence, which is what the calibration was designed to refuse. Gorilla chr20 shows the same pattern (181 of 213
stopped by the calibration). ⚠ A relaxed rule (for example best-hit family + counting hit + the 0.60 floor, without
the calibration) is a different rule: it needs its own pre-registration, development on chr16, and a verdict on
untouched data (chr6 is now used).

**Prediction scorecard.** Predicted ⚠ SAFE, SMALL with sensitivity < 0.10 because most reference pairs would lack a
protein hit (r1026/r1101). The verdict matched, the mechanism did not: 27 of 30 missing members DO have a counting
protein hit to the right family; the binding clause is the pre-registered calibration, not the search.

**Not yet measured.** Liftoff (no species merged yet), the genome-wide per-sample rows (the driver stage
`families-protein` has not run; human_testis repeats the Compara verdict with chr16 dropped).
