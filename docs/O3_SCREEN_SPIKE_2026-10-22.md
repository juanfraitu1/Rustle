EXPLORATORY (spike, post hoc)

# O3 RNA-only "look harder" screen, tiers 1-2: a 45-minute spike on HG002 and HG00733 (2026-10-22)

**What.** Tiers 1 and 2 of `docs/O3_FINGERPRINT_SURVEY_2026-10-22.md` (definitions: `docs/o3survey_2026-10-22/C_certificate.md` §1-3).
- Built as `bench/o3_graph/o3screen.py OUTDIR COPIES GENOME OUT.tsv`, one row per family (2 sanity tests: `test_o3screen.py`). It
  reads only the frozen O2-VC runs `/mnt/linuxdisk/tmp/o3_graph/o2vc/{hg002,hg00733}`: HG002 about 2.3 min and 220 MB, HG00733 4 s.
- The truth is used only for the family labels. Every variant is reported, including those made after looking (marked †).

## Definitions as built (constants: a fact seen by >= 2 reads, a conflict at >= 2 positions)
- **Informative positions of copy c:** every position of c in the family's `columns.tsv`, any kind, any tied read.
- **Pattern:** a read's bases at the informative positions it covers, base quality >= 10, in c's genome strand.
  - Teachers come from `teach/<family>.bam` (attributed to copies as `o2vc_data.scan` does); tied reads come from `columns.tsv`
    (the strand is told from the genome base against the row's `ref`).
  - A read covering no position has no pattern.
- **Conflict:** two patterns differ at >= 2 positions that both cover.
- **Classes: greedy first fit,** not exact match.
  - Distinct patterns are taken most positions first, then most reads; each joins the first class with no member it conflicts with.
  - A class of >= 2 reads is supported.
  - A class's pattern is its certified consensus: positions where >= 2 members carry the same base and it is the strict majority.
    The brief left this open.
- **Tier 1 (per copy):** a greedy multi-start clique of pairwise-conflicting supported teacher classes. It is a lower bound on
  the clique number. A copy is flagged when the clique is larger than its capacity.
  - Capacity: 1 on a non-PAR chrX/chrY copy of a male (PAR from `hg002truth.PAR`), else 2.
  - Sex was checked in the data. HG002 is male: tied reads in chrY's HSFY1; 25 of its 54 chrX/Y copies get capacity 1. HG00733
    is female: 0 of 127 non-PAR chrY copies have an ortholog, and 28 non-PAR chrX copies have both a mat and a pat one.
- **Tier 2 (per family):** the tied reads of `calls.tsv`, classed per anchor copy and counted once per read set.
  - A supported tied class is **novel** when it conflicts with every supported teacher class of every copy its reads tie to; each
    copy is compared in its own coordinates.
  - **t2 (literal):** a copy with no supported teacher class imposes nothing (vacuous).
  - **t2t (taught):** every candidate copy must have >= 1 supported teacher class.
- **Today's flag:** >= 1 row in `<family>.flags.tsv`.

**Variants** (flags of the script):
- **v1** is the default above.
- `--exclude-editing` drops `editing` rows, as C §1 drops A>G.
- **v2 `--pair-certified`:** the two conflicting positions must be carried *together* by >= 2 members of each class (C S1).
- **v3 `--merge-single`†:** classes one position apart are merged by single linkage (C §1). Added after diagnosing v1's cliques.
- **v4 = v3 + `--leak-closure`†:** tier 1 needs a Hall deficiency > 0 (C S4). Each clique class may also take the slots of any
  other copy with a teacher class it does not conflict with. Added after the HG00733 false alarms were traced to leaked teachers.
  - The comparison is made in the other copy's coordinates, through the `columns.tsv` homology.
  - A copy that no tied read links to c shares no position, so it counts as compatible (the conservative choice).
    `--leak-needs-overlap` (v4o†) requires >= 1 shared position instead.
  - A first leak-closure run compared positions without translating them, so it was vacuous. It was fixed before every number here.

## HG002 (94 families with tied reads; the 15 PAR families excluded; 73 "others")
Positives: (a) SPDYE1, RASA4B, USP17L11, POTEM, BMS1P1; (b) = (a) + NPIPB2.

| arm | (a) | (b) | others flagged (false alarms) |
|---|---|---|---|
| tier 1 v1 | 2/5 SPDYE1, RASA4B | 3/6 +NPIPB2 | 9/73: AC1275021, ANKRD36C, DUSP22, GTF2H2C, NBPF1, NOMO1, SMG1, SMN2, THOC3 |
| tier 1 v1, no editing | 1/5 SPDYE1 | 2/6 | 8/73 (SMG1 drops) |
| tier 1 v2 | 2/5 | 3/6 | 9/73 (same as v1) |
| tier 1 v3† | 1/5 SPDYE1 (3/2) | 2/6 +NPIPB2 | 5/73: AC1275021, ANKRD36C, GTF2H2C, NBPF1, SMG1 |
| tier 1 v1 + leak† (overlap required) | 1/5 SPDYE1 (1/5) | 1/6 (2/6 +NPIPB2) | 6/73: DUSP22, GTF2H2C, NOMO1, SMG1, SMN2, THOC3 (7/73 +NBPF1) |
| tier 1 v4† | 0/5 | 0/6 | 2/73: GTF2H2C (.1 deficiency 3), SMG1 (.2, 1) |
| tier 1 v4o† | 0/5 | 1/6 NPIPB2 (.9, deficiency 1) | 3/73: GTF2H2C, NBPF1, SMG1 |
| tier 2 v1 literal / taught | 0/5 / 0/5 | 1/6 NPIPB2 (31 reads) / 0/6 | 5/73 HSPA6, MSTRG3150, SMN2, THOC3, TRIM74 / 4/73 |
| tier 2 v2 literal / taught | 0/5 / 0/5 | 1/6 / 0/6 | 2/73 MSTRG3150, SMN2 / 1/73 SMN2 |
| tier 2 v3† (= v4) literal / taught | 1/5 RASA4B (119 reads) / 1/5 | 2/6 +NPIPB2 / 1/6 | 7/73 DUSP22, GOLGA8A, HSPA6, MSTRG3150, SMN2, THOC3, TRIM74 / 6/73 (no MSTRG3150) |
| today's flag | 1/5 SPDYE1 | 2/6 +NPIPB2 | 7/73: AC1275021, GOLGA8A, MSTRG3150, NCF1, NOMO1, SMG1, SMN2 |

The PAR families flagged (excluded from the table) were CD99 and SLC25A6 (tier 2) and SLC25A6 (today's flag).

**Ranking by tier-2 read support (v1, literal).**
- THOC3 2,887 > SMN2 1,594 > MSTRG3150 144 > HSPA6 109 > **NPIPB2 31** > TRIM74 20. Every other family scores 0.
- AUC-like score (positives against the 73 others, ties count 0.5):
  - tier-2 reads: 0.47 (a) / 0.55 (b);
  - today's flag reads: 0.56 (a) / 0.64 (b).
  - With nearly every family at 0, both are close to chance.
- **v3 ranking:** THOC3 2,892 > SMN2 1,602 > DUSP22 152 > GOLGA8A 146 > MSTRG3150 144 > **RASA4B 119** > HSPA6 109 > **NPIPB2 31** >
  TRIM74 20, then 0.
  - AUC-like score: tier-2 reads 0.55 (a) / 0.61 (b); today's flag 0.56 / 0.64.
  - SPDYE1 is caught only by tier 1 and today's flag.

**Why three positives are missed in every version.**
- USP17L11 and POTEM have **no taught copy**: 0 copies with a supported teacher class, from 12 and 26 teachers over 62 and 11 copies.
- BMS1P1 has 2 taught copies, 13 tied reads, and neither a clique > 2 nor a recurrent tied class.

## HG00733 (simulated, female, no donor-only copies: every flag is a false alarm; 174 families with tied reads)
The donor lacks a reference copy in 46 families, reproduced from `orthologs.tsv` against `pool_copies`.

| arm | false alarms / 174 | in the 46 absent-copy families |
|---|---|---|
| tier 1 v1, v1 no editing, v2, v3 | 3: GOLGA6L6, PRAMEF10, SPATA31A1 | 1/46 SPATA31A1 |
| tier 1 v1 + leak†, v4† | **0** | 0/46 |
| tier 1 v1 + leak†, v4o† (overlap required) | 1: GOLGA6L6 (.1, deficiency 1) | 0/46 |
| tier 2 v1 literal / taught | 10: FFAR3, KRT17P1, LIMS1, LYZL1, NCF1, OR1E1, PRAMEF11, RHOXF2B, TPSB2, ZC3H11A / 8 | 2/46 LIMS1, PRAMEF11 / 2/46 |
| tier 2 v1 no editing, literal / taught | 8 / 6 | 1/46 PRAMEF11 |
| tier 2 v2 literal / taught | 5: FFAR3, OR1E1, PRAMEF11, TPSB2, ZC3H11A / 5 | 1/46 PRAMEF11 |
| tier 2 v3† (= v4) literal / taught | 8 (v1's minus PRAMEF11, ZC3H11A) / 6 | 1/46 LIMS1 |
| today's flag | 5: FFAR3, GOLGA6L6, KRTAP21, TPSB2, ZC3H11A | 0/46 |

## Diagnosis (post hoc; it explains the variants)
- **v1's cliques are mostly fragmentation.**
  - THOC3.0 has 4,054 teachers and 508 informative positions. They give 109 supported classes: one of 3,685 reads and 108 of 2-9
    reads. The greedy clique is 33.
  - Of the 288 reads in the small classes, 267 differ from the main consensus at exactly 1 position.
  - So a large class sheds reads, because "non-conflicting with every member" adds two members' independent errors together.
    Each fragment then needs only one certified deviation to conflict with every other fragment.
  - v2 does not stop this, since each fragment carries both of its positions; v3 does.
- **HG00733's tier-1 false alarms are teacher leakage** (traced with `sim/truth.tsv`). The third class of every flagged copy is reads
  from a paralog that minimap2 placed as teachers:
  - GOLGA6L6.1: 2 reads from .3;
  - PRAMEF10.1 and .4: 5 and 2 reads from each other;
  - SPATA31A1.1 and .3: 2 and 4 reads from each other.
  Leak closure removes all three. GOLGA6L6 goes only through the vacuous rule: no tied read links .1 to .3. Requiring a shared
  position leaves it (v4o: 1/174).
- **Leak closure is per copy here.** Another copy's slots are not reserved for that copy's own haplotypes, so it is weaker than
  C's family-wide Hall test.
  - SPDYE1.18 is lost in v4 and in v4o: at least one of its 3 classes shares positions with a sibling copy's class without
    conflicting.
  - NPIPB2.9 survives only in v4o.
  - The family-wide matching, with the same haplotype seen from two copies counted once, was not built.
- **RASA4B.2** (v1) is two real haplotypes plus a 13-read class at two adjacent `editing` positions.

## Reading
1. As specified (v1), neither tier beats today's flag. Tier 1 catches 2/5 but false-alarms in 9/73 HG002 and 3/174 HG00733 families (bar: 0); tier 2 catches 0/5. Its cliques are fragmentation and teacher leakage, not haplotypes.
2. C's own rules, added post hoc† (merge, then leak closure), make tier 1 specific (v4: 0/174 HG00733, 2/73 HG002) but blind (0/5; v4o 1/6 via NPIPB2). Tier 2 v3 catches RASA4B (and NPIPB2 in (b)) but ranks CNV-prone SMN2, GOLGA8A and THOC3 above them (AUC-like 0.55 / 0.61 vs today's flag 0.56 / 0.64).
3. USP17L11, POTEM and BMS1P1 have no taught copy or too few reads, so no teacher-based tier can see them. Nothing here is claimable; next is a family-wide Hall test and a prereg on a fresh donor with simulated donor-only copies. Absent copies stay invisible (0/46).
