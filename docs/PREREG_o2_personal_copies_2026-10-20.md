# Pre-registration: personal copies against the corrector, on gorilla KB3781's own reads with exact truth (2026-10-20)

**User (2026-10-20):**
- "lets go back to O2, and assume we have a family with members with different degrees of similarity (90 to slighly less
  than 100%); what is still missing to be able to disambibuate (ambiguous AS) to know the real provenance and use that for
  assembly"
- On the proposal to test per-copy sample sequences against the corrector's pileup: "yes please, also maybe we can use
  gorilla for this since we have the fibroblasts and T2T genome?"
- On the design below: "yes please proceed".

## The question

The corrector (`o2correct.py donor`, registers 1343-1346) learns each copy's sample alleles as a **pileup**: per copy, how
many of its confident reads carry each substitution or indel. A tied read's events are then checked one by one against
that pileup.

A **personal copy** is the copy's own sequence in the sample: the reference copy with the events its confident reads
carry written in, split into at most two haplotypes by the reads that carry them. A tied read is realigned to the personal
copies of its tied copies.

The personal copy keeps what the pileup loses: which events sit together on one haplotype, and the read's alignment to the
sample's sequence rather than to the reference.

The experiment asks:
1. Do personal copies place tied reads better than the pileup?
2. How close does either get to the animal's real copies?

## Substrate: one animal, his reads and his own genome

- **The reads.** KB3781 ("Jim", the animal of the mGorGor1 T2T assembly): fibroblast Iso-Seq, SRR27178662, 13,516,443
  molecules and 34,762,361 records. They were aligned to the RefSeq primary GCF_029281585.2 (`GGO.fasta`) by minimap2 2.31
  `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`
  (`/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam`).
- **The reference is a mosaic of his two haplotypes** (`rna_allele/chrmap.tsv`, sha1 e9182791, each chromosome checked
  identical in sequence):
  - his father's: chromosomes 1, 2 (NC_086017.1), 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21 and Y;
  - his mother's: 5, 6, 10, 11, 17, 20, 22, 23 (NC_086018.1) and X.
  - This corrects the design note in chat, which said the autosomes were all his father's.
- **The ALT haplotype** is the other parent's copy of each autosome, plus every unplaced contig of either assembly. His
  assemblies are GCA_028885475.2 (paternal) and GCA_028885495.2 (maternal), in `gorilla_haps/{pat,mat}.fa`.
- **Two built-in groups of reads.** A read from the haplotype the reference was taken from has a sample copy identical to
  the reference: the **one-genome control**. A read from the ALT haplotype carries the other parent's alleles: the **test**.

## The families

- **Catalog.** The KB3781 fibroblast catalog GWFAM: `rna_allele/fibro_cat.copies.tsv`, sha1 7521851f, 666 families, 3,981
  copies.
- **Rule.** A family qualifies when every pair of its members' main transcripts aligns at identity in [0.90, 1.00), covering
  at least 50% of the shorter. The test is `minimap2 -X -c -x asm20` with `o2perfamily.qualifies`, the human list's rule.
  - It is run one family at a time (`o2personal.py ava`), since only same-family pairs count.
  - A family of more than 20 members is first tested on 20 of them; a pair failing there fails the family.
- **Result.** 283 families qualify, with 716 copies (`personal/qual.json`, sha1 f4bab888; `qual.copies.tsv`, sha1 78a3afa8).
  NPIP does not qualify.
- **The hold-out is by family.** The parity of sha1(family name), fixed before any result was seen (`o2personal.half`):
  - development: 140 families;
  - held-out: 143 families.

  | minimum pair identity | development | held-out |
  |---|---|---|
  | [0.90, 0.95) | 15 | 11 |
  | [0.95, 0.98) | 50 | 64 |
  | [0.98, 0.99) | 37 | 33 |
  | [0.99, 1.00) | 38 | 35 |

## The domain (the corrector's, unchanged)

A read is in the domain when:
- its alignments within 10 AS of its best lie on two or more copies of one qualifying family;
- none of them lies on a locus off those copies;
- minimap2's primary is among them.

The reads are every alignment of the 450,984 reads with any alignment on a qualifying copy (`personal/sub.bam`, 4,055,491
records, `names.txt` sha1 26dd6d3c).

## The truth: the read's best placement in the animal's own genome

**Alignment.** Each domain read is aligned to his maternal and his paternal assembly (`mGorGor1.{mat,pat}.splice.mmi`):
minimap2 2.30, `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

**Each alignment maps to a reference copy:**
- On a haplotype contig the reference is identical to, it maps to the copy it overlaps there, in the reference's
  coordinates (`summary.attribute`, the corrector's own attribution).
- On an ALT contig, it maps to the copy whose ALT ortholog it overlaps. A copy's ALT ortholog:
  - its reference span ± 20 kb, aligned with `minimap2 -x asm20 -c` to the other haplotype's whole assembly;
  - the longest primary alignment with MAPQ >= 20 through which its exons lift (`donorsim.pick_ortholog`);
  - on the homologous chromosome or an unplaced contig;
  - two copies whose orthologs overlap are both dropped (`donorsim.drop_ambiguous`).

**The truth** (`o2personal.truth_of`) is the one copy that every alignment within 5 AS of the read's best lies on. 5 AS is
one substitution under splice:hq (match 1 + mismatch 4).

| status | meaning |
|---|---|
| determined, **refcons** | The reference's own copy is among them: the sample's copy equals the reference over the read. |
| determined, **altonly** | Only the ALT haplotype's copy: the read carries the other parent's alleles. |
| undetermined | Two copies, or a copy and an off-copy locus, lie within 5 AS. No method can place the read. |
| off-copy | The best placement is on no reference copy: a copy the reference lacks, an unannotated locus, or an ortholog that did not lift. |

**Scored reads (S):** determined, with the truth among the read's tied copies.

## The arms

Each arm gives one copy per read. An arm that abstains keeps minimap2's primary.

1. **minimap2**: its primary.
2. **Corrector**: `o2correct.py donor`, version 2 without `--guard`, unchanged since bfa15f12.
3. **Personal copies, phased** (`o2personal.py`). Per copy:
   - **Reads:** the corrector's own confident reads (primary on the copy, every alignment on another copy of its family more
     than 9 AS lower).
   - **Events:** those the corrector calls explained (>= 2 reads, >= 10% of a depth >= 3). Events at >= 80% of the depth are
     **carried**.
   - **Haplotypes:** the other events split into at most two haplotypes by a 2-means over the reads (`o2personal.phase`).
   - **Sequence:** each haplotype is the reference span ± 5 kb with its events written in.
   - **Decision:** each tied read is realigned to the personal copies (the same minimap2 options, `-N 200`). The tied copy
     with the unique best AS takes it. A shared best, or a tied copy the read does not align to, keeps minimap2's pick.
4. **Oracle**: the same decision over the copies' true sequences: the reference span ± 5 kb and its ALT ortholog ± 5 kb.
5. **Exploratory, no bar**: personal copies unphased, one sequence per copy (the events more than half the depth carries).

## Scoring (`o2personal.py score`)

**Per arm, on S:**
- right and accuracy;
- against minimap2: **fixes** (minimap2 wrong, arm right) and **breaks** (minimap2 right, arm wrong);
- reads decided, and decided right.

**Strata:**
- refcons and altonly;
- exact ties (AS margin 0 between the two best tied copies) and near ties;
- the family's minimum pair identity, in the four bins above;
- each family.

**Also reported:** undetermined, off-copy, and truth outside the tie set.

## Bars, held-out half

- **G0 (gate):** oracle accuracy on S >= 0.99. Failing it means a pipeline error: INVALID.
- **P (power):** |S| >= 300, from >= 11 families, with minimap2 wrong on >= 50 of them.

**V1. Personal copies beat the corrector:**
- **G1:** personal accuracy >= the corrector's + 1.0 point.
- **S1:** personal fixes >= 5 × personal breaks.
- **F1:** >= 11 families where personal is right more often than the corrector, and no more than a fifth as many (rounded
  up) where it is right less often.
- **H1:** on refcons reads, personal accuracy >= minimap2's − 0.2 points.

**V2. The corrector beats minimap2** (secondary; its first real-read test with exact truth):
- **G2:** corrector accuracy >= minimap2's + 2.0 points.
- **S2:** fixes >= 5 × breaks.
- **F2:** >= 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up) where it breaks
  more than it fixes.
- **H2:** on refcons reads, corrector accuracy >= minimap2's − 0.2 points.

## Verdict

- **INVALID:** G0 fails.
- **UNDERPOWERED** (both V1 and V2): G0 holds and P fails.
- **PASS:** all four of that verdict's bars hold. **FAIL:** any of them does not.

**A V1 PASS would claim:** on one gorilla's own fibroblast reads, with exact truth from his own assemblies, the sample's
phased copy sequences place tied reads better than the corrector's pileup. This holds in more than 10 near-identical
families and leaves reads whose copy equals the reference unharmed.

**It would not claim:**
- other tissues or animals;
- copies absent from the reference (off-copy: O3);
- reads tied in the sample itself (undetermined);
- a better assembly (no second round here).

## Development, then the freeze

The development half (140 families) may be run and looked at freely. On it, I may:
- fix bugs;
- change the personal copies' construction, the realignment's minimap2 options and FLANK_P;
- change the truth's mechanics (contig names, ortholog lift).

Fixed for the whole experiment:
- the truth margin (5 AS);
- the domain;
- the corrector (bfa15f12);
- the bars.

Every change goes into Amendment 1, with its reason, at a commit. The held-out half is scored once, with the code at that
commit.

## Rules

- Gorilla numbers are never pooled with human ones.
- Heavy jobs run in the foreground, one at a time: the truth alignments (about 14 GB of index each) and the ALT ortholog
  alignments.
- The metric-auditor checks the held-out numbers before any verdict is reported.
- Not run, and kept as a later replication: the mirror, with reference = the other haplotype.

## Amendments

(none yet)
