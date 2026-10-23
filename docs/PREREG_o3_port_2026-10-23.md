# Pre-registration: O3 port (Strainy's phasing + IsoCon's test), gorilla KB3781 (2026-10-23)

Committed before the port has been run on any panel. Design: `docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md`.
The bars are the shipped stage's own bars (`candidate_copies`, Amendments 13-15: register rows 1221-1231 and 1247), so that the
port can replace the stage on equal terms.

## Arms, the same reads and the same scorer
- **P:** the port, run as designed.
- **K:** `candidate_copies` as it ran under Amendment 15. Its products are already on disk: a13a15 / a14a15 / a15h. It is not
  re-run.

## Decision rules
- **P1, development (deletion panel, a13a15).**
  - **Pass:** D right >= 0.80 x C = 5,842 (C = 7,303, Amendment 13b's comparator) AND false moves <= 5% of S reads.
  - **K on the same panel:** 7,882 and 0.77%, which passes.
- **P2, development (no-deletion control, a14a15).**
  - **Pass:** families with >= 1 false flag (class b / c / pri against mat / pat, Amendment 9's rule) <= 8 of 53.
  - **K:** 10 of 53, which fails.
- **P3, held out (a15h; run ONCE, after P1 and P2, with frozen code).**
  - **Pass:** families with >= 1 false flag <= 4 of 30 (H1) AND false moves <= 5% (H3).
  - **K:** 6 of 30 and 0.26%; H1 fails.
- **Adopt the port as the O3 stage iff P1, P2 and P3 all pass.**
  - P passes P1 but fails P2: "detects but not specific", not adopted.
  - P fails P1: "not sensitive", not adopted.

## Rules for development
- P1 and P2 are development. A change after looking is allowed only if it adds no new constant, is recorded in the doc with its
  reason, and is re-run on BOTH P1 and P2.
- Then the code is frozen (commit hash written here) and P3 is run once.

## Reported, not decided
- The per-copy haplotype count in the control against mat / pat: allele-level truth for every copy.
- How many certified flags depend on an editing-type column.
- Runtime.
- HG002 SPDYE1 (development case, human).

## Development log (P1 / P2 are development; every change below adds no new constant)
- **v1** (1cb54fa7), on GWFAM175 / GWFAM4 only:
  - per-read empirical error rates (NM / aligned bases, IsoCon) instead of PacBio QVs, which are capped at Q93 and let 2 carriers
    of a recurrent error pass;
  - a conflict needs >= 2 differing columns.
- **v1 deletion panel** (53 families), descriptive: 91 certified haplotypes in 35 families. Of their clusters, 40 hold only D
  reads, 47 only S reads and 4 both. 16 of the 53 deleted copies are caught: in 31 of the 53, no D read aligns to a survivor
  (A11).
- **v1.1:** conflicts are proven only on pile-level columns. Columns that pass only inside a refined cluster can come from
  recurrent errors or chimeras.
- **v1.2:** an extra copy needs >= 2 private alleles relative to the slot holders (mosaics excluded).
- **P2 on v1.2** (full classification, `ACC=a14 A=p14`): 53 flags in 22 of 53 families.
  - Classes: c_unmatched 45 (mostly 2-7 reads), b_allele 5 (219 / 118 / 97 / 37 / 3 reads), pri 2, a_haplotype_only 1
    (GWFAM175, 181 reads).
  - **C1' 22/53, FAILS** (bar <= 8). The shipped stage scored 10/53.
- **Diagnosis.** Of the surplus haplotypes, 4 of 72 are made of short reads (median read length 0.79-0.85 of the pile's p90),
  so requiring full-length reads would not remove them. Their private alleles are transition-rich and none are all A>G, so they
  are not editing. 28 of the 72 are mosaics of the copy's two alleles.
- **v2:** the family-wide DADA2 rule runs before the slots. A haplotype that equals a more abundant haplotype, or is a left/right
  join of two parents, takes no slot. Its possible parents are every copy's reference and the more abundant haplotypes, mapped
  by asm20 homology. Result: P2 19/53, P1 13 families with D.
- **v3:** MIN_SUPPORT 6, the shipped stage's own flag floor (candidate_copies --min-support 6, ruling R1), so the comparison is
  on equal terms: a haplotype needs 6 reads to take a slot or be certified.
  - Slot order compared on development: by divergence, P1 has 14 clusters, 1,031 D reads and 2 S reads; by read support, 301 D
    and 396 S. Divergence was kept.
  - P2 had 0 flags under both orders. GWFAM175's real paternal-only copy is no longer certified: its reads tie to a copy with a
    free slot, which is Theorem A's ambiguity.
- **v4 = FROZEN (1b4535cc):** v3 plus the IsoCon side, i.e. the stage's flagged candidates whose members are >= half via-B (unaligned)
  reads. On P1 that adds 13 candidates in 12 families (804 D, 5 S reads); on P2, 0.
- **P1 (v4, harness ACC=a13 A=p13): D right 3,995 (54.7% of C) < 5,842: FAILS.** False moves 7 / 41,727 = 0.02% (K: 7,882 and
  0.77%).
- **P2 (v4): 0 flags in 0 / 53 families: PASSES** (K: 10 / 53).
- **Verdict so far: "detects but not sensitive", not adopted** (pre-registered reading of a P1 failure). P3 is run once on the
  frozen code to record held-out specificity (H1 / H3); it cannot change the verdict.
- **P3 (frozen 1b4535cc, run once, harness ACC=a15h A=p15h):** 6 flags in 4 / 30 families.
  - **H1 3 / 30, PASSES** (<= 4). Classes: 1 a_haplotype_only (GWFAM390, 8 reads, maternal-only locus, identity 1.0000),
    3 b_allele, 2 c_unmatched.
  - **H3 178 / 27,757 = 0.64%, HOLDS** (<= 5%). 362 reads go to their own candidate.
- **FINAL: P1 fails, P2 and P3 pass. Not adopted ("detects but not sensitive").** Register row added.
