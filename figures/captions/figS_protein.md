# Supplementary Figure S-P: the manual extra-sensitive step (protein attachment)

**Status.** Pre-registered in `docs/archive/2026-09/PREREG_protein_attach_2026-09-25.md` before any number of the step existed.
Development-scope tables (PROVISIONAL): the de novo families of the Fig. 7 runs on four chromosomes. Genome-wide rows
appear per sample once the driver's `families-protein` stage has written `${work}/runs/<id>/<id>.fam_protein.*`;
Liftoff rows appear once the Fig. 8 self-lift of a species is merged. `python3 figures/fig_protein_supp.py data`
rebuilds the tables (the development runs are light: 1-26 s each), `plot` renders, `summary` prints every number
quoted here.

**Claim.** Protein is not part of the default de novo family definition. Run by hand after the default families
exist, the step attaches very few loci: 0 of 1,470 candidates on human chr16 (development), 1 of 2,339 on human chr6
(held out, never used), 1 of 1,007 on gorilla chr20 (NC_073244.2) and 0 of 887 on gorilla chr10 (NC_073234.2); the
shuffled-protein null attaches none anywhere. The chr6 attachment is correct against Ensembl Compara (paralogues of
any age; 1 of 1), as are all 14 of chr6's judgeable family members, but it attaches 1 of the 30 missing members
Compara identifies. Pre-registered verdict: ⚠ SAFE, SMALL (as predicted), on a single judgeable attachment.

## Caption

**Supplementary Figure S-P | What the optional protein step adds to Rustle's default de novo families, scored against
external references.** The default de novo families are unchanged: IsoSeq reads are assembled into loci (seeded with
secondary alignments within 2% of the molecule's best alignment score); each locus is represented by the exons of its
most-read transcript; loci are grouped when their genomic sequences align (minimap2, ≥ 70% identity) and share ≥ 60%
of the smaller locus's exonic sequence; Markov clustering (inflation 2.8). The protein step is run by hand
(`tools/protein_attach.py`; driver stage `families-protein`, never part of `all`). It never changes a default family
and never merges two families; it can only add a locus the default rule left out to one existing family.

*The step.* Each locus's representative exons are spliced from the genome and translated: the longest stop-to-stop
reading frame on the transcribed strand (no start codon required; an `N` codon ends a frame), among frames less than
half soft-masked (a frame made mostly of transposable-element sequence is skipped), kept if ≥ 100 amino acids. Every
locus protein is searched with BLASTP 2.17.0 (E ≤ 1e-5) against the proteins of the family members. A hit counts when
its non-overlapping alignments, chosen by bit score, cover ≥ 30% of the longer protein; its identity is taken over
those alignments. Each family is calibrated by its own members: its loosest member's best identity to another member
(a family whose members share no protein homology receives nothing). A locus the default rule left out joins the
family of its best BLASTP hit only if it is at least as identical to that family as the family's loosest member,
at least 60% identical (the protein identity below which RNA of paralogues no longer aligns), and no other family
passes the same test. A locus whose exons overlap a member's is never a candidate. Null: every candidate protein
shuffled (fixed seed), same search and rule.

**a**, Loci attached per substrate (bar) and null attachments (open circle); the text gives candidates, families
receiving loci, and pairs of families whose members pass the same test (reported as merge proposals, never applied).
**b**, Candidates with a BLASTP hit to a family member, by the first clause that stopped them. **c**, Precision
against each reference: open circles, the default families' own members (a member is correct when the reference
relates its gene to the gene of another member of its family); filled circles, attached loci (correct when the
reference relates its gene to a member's gene); lines, 95% Wilson intervals; right, correct / judgeable. **d**,
Missing members (human): loci the default rule left out whose gene no family holds and which Compara relates to a
member's gene; how many have a protein, a counting hit to that family, and were attached to it.

*References.* Ensembl Compara release 116 (human only): two genes are related when a chain of Compara paralogue pairs
joins them; "any age" uses every duplication node (one gene tree; separates paralogues from proteins that merely share
a fold), "primates" only duplications within primates. Liftoff copies (Fig. 8 self-lift, `-copies -sc 0.95`; a locus
matches a placement covering ≥ 50% of its exonic bases) can only certify an attachment, since Liftoff lifts copies
≥ 95% identical; not yet available. *Protein-homology families (grey, asterisk) are built from protein homology of
the annotation, so agreement with them is partly by construction; they are a secondary reference only. A locus takes
the gene (RefSeq gene or pseudogene record) whose exons cover ≥ 50% of its exonic bases; otherwise it is not
judgeable. Species are never pooled.

*Substrates and exposure.* Human A119b (T2T-CHM13 v2.0) chr16: development (family rules were tuned here; the step's
constants were fixed before it was run). Human chr6: held out, never used for a family or protein decision (scored
once for Fig. 7). Gorilla OR6737, testis (mGorGor1, GCF_029281585.2) chr20 (NC_073244.2): the chromosome the seeding
rule was chosen on; chr10 (NC_073234.2): held out, never used. Each is the genome-wide assembly restricted to that
chromosome, grouped by the default family rule (the Fig. 7 development runs).

## Numbers (from `figS_protein_substrates`, `figS_protein_refs`)

| substrate | loci | families (calibrated) | loci with a protein | candidates | with a member hit | attached | null | merge proposals |
|---|---|---|---|---|---|---|---|---|
| human chr16 (development) | 2,802 | 111 (51) | 1,746 | 1,470 | 28 | 0 | 0 | 2 |
| human chr6 (held out) | 5,012 | 54 (11) | 2,456 | 2,339 | 84 | 1 | 0 | 0 |
| gorilla chr20 (NC_073244.2) | 1,142 | 29 (24) | 1,087 | 1,007 | 213 | 1 | 0 | 15 |
| gorilla chr10 (NC_073234.2) | 956 | 6 (1) | 895 | 887 | 1 | 0 | 0 | 0 |

- Why candidates with a member hit were not attached (panel b): chr16 9 less identical than the family's loosest
  member / 9 hit < 30% of the longer protein / 10 family without protein homology among its members; chr6 29 / 6 / 48
  (+ 1 attached); gorilla chr20 181 / 8 / 10, 11 below 60% identity, 2 passing for two families; gorilla chr10 1
  uncalibrated family.
- chr6 attachment: the locus of HLA-E, to the HLA class I family (identity 0.766, 89% of the longer protein; the
  family's loosest member 0.750); a sub-threshold nucleotide alignment to the family also exists. Correct against
  Compara at any age; not a duplication within primates (0 of 1 at that level).
- Default members, Compara any age: chr16 31 of 31, chr6 14 of 14 (95% CI 0.78-1.00); within primates chr16 31 of
  31, chr6 8 of 14.
- Missing members, Compara any age: chr16 8 (8 with a protein, 4 with a counting hit to the related family, 0
  attached); chr6 30 (29, 27, 1): sensitivity 0.033. Within primates: none on either chromosome.
- chr6's 27 missing members with a counting hit are all hit best in the family Compara relates them to; 20 fall below
  that family's loosest-member identity and 6 hit a family whose members share no protein homology (HLA class II,
  MICA/MICB, ULBP1-3, BTN2A1/2A2/3A1, GSTA3/4, TUBB/TUBE1, SLC22A2/3/16, ZNF311/391). The family calibration is the
  clause that binds (observed after the verdict; changing it needs its own pre-registration and a substrate other than
  chr6).
- Merge proposals (never applied): chr16 2 (Compara: 0 judgeable); gorilla chr20 15 (protein-homology families relate
  all 15, circular).
- Proteins: 442 of chr16's loci (1,161 on chr6) have their only ≥ 100 aa frames mostly soft-masked; masking differs
  between genomes (gorilla chr20: 37).

## Caveats

- One judgeable attachment on the held-out chromosome: the verdict is a scope statement, not a precision estimate.
- Liftoff (the only external reference for the apes) is not merged yet; gorilla rows have only the circular
  secondary reference.
- The development tables are single chromosomes of the genome-wide assemblies; the per-sample genome-wide rows
  (human testis, gorilla KB3781, chimpanzee PTR, orangutan PPY held out) replace them when the driver stage has run.
- Compara judges only protein-coding genes present in its paralogue table under the same symbol; pseudogene copies,
  the step's other target, are judgeable only through Liftoff.
