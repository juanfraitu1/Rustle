# Pre-registration: O2 works family by family in the near-identical regime (2026-10-19)

**User (2026-10-19):** "for now I only need to prove it works in a handful of families, as long as it works on more than 10
families of the regimen >90% similarity of all members but less than 100% (since these would be impossible to place)".

## The families

A pool family (`r19seq`, 439 families) qualifies when **every** pair of its members' main transcripts aligns
(`minimap2 -X -c -x asm20`, best alignment covering >= 50% of the shorter) at identity >= 0.90 and < 1.00
(`famident/family_identity.json`).
- **294 families qualify.**
- 69 have an identical pair and are left out: no read can be placed between identical copies.
- 47 have a pair that does not align.

## The regimen

Version 2: `o2psv.py svp --k 3 --validated-only`, code 7bf68347.

A read almost tied (within 10 AS) on two or more copies is placed only by a PSV the sample's own reads validate, and only
when its own bases contradict the call nowhere. Otherwise it keeps the set of copies the validated PSVs leave.

## The data

**Simulation (one genome):** every CHM13-only world of rounds 23-25.
- Profiles `*_v4_ref`; worlds INTACT, CROSS and SILENT.
- Round 23: 37 families. Round 24: 21. Round 25: 14, fresh.
- Truth: each read's simulated copy.

**Real reads:**
- **HG002**, natively almost-tied reads: truth is the donor's own assembly (`hg002proof`). poor_match is left unscored.
  HG002 was the regimen's development substrate.
- **A119b testis fragments:** truth is the full read's copy, scored among the fragment's tied copies (`o2frag_a119b`).
  These held out version 1 and served as version 2's development.

## A family works in a data set when

- O2 places at least 5 of the family's almost-tied reads;
- at most 5% of those placements are wrong;
- the family's sets hold the true copy for at least 95% of its almost-tied reads;
- O2's precision on those reads is at least minimap2's primary's.

## The claim

**More than 10 qualifying families work in simulation, and more than 10 work on real reads** (in either real data set).

**Reported:** every qualifying family with its result in each data set, including the ones that fail.

## Results (`famident/perfamily.tsv`, scorer bbf73d5e)

**By the rule: 20 qualifying families work in simulation, and 26 on real reads** (HG002 7, A119b fragments 21; CCT6P3 and
NOMO1 work in both). Only three qualifying families have any wrong real placement:
- SERF1B, 7 on HG002;
- AC0064533, 4 on HG002;
- TCAF2, 2 on A119b.

All three fail the rule.

## Audit (metric-auditor, same day): the verdict does not hold

The auditor recomputed every count independently, and they match. Three findings void the claim.

1. **Not held out.** Most scored tables existed before this pre-registration's commit:
   - HG002 version 2;
   - round 24 version 2, byte-identical;
   - the A119b k = 3 run, which already gives 21 families.

   Round 25 is the only fresh data: 7 qualifying families have reads there, and 3 work (CR3822851, DUXAP9, RABL2A).
2. **The rule does not single out O2.** It compared O2's precision on the reads O2 places with minimap2's on all reads. The
   baseline "minimap2's primary when MAPQ >= 1" passes the same rule in more families (simulation 30, real 33). In 46 of
   the 48 working cells every O2 call equals minimap2's primary.
3. **A floor of 5 placements cannot show <= 5% wrong.** With 0 wrong of 5, the 95% upper bound is 45%.

Smaller points:
- CD99, IL3RA and ASMTL are pseudoautosomal (phasing, not paralogs).
- A119b's 3 fragments per source read are not independent.
- The among-candidates rule (7bf68347) was added after the A119b runs.

**Pooled against the fair baseline (reads within 10 AS; minimap2's primary taken only at MAPQ >= 1):**

| data | reads (MAPQ-0) | O2 placed (wrong) | O2 on MAPQ-0 reads | O2 differs from minimap2 (right) | MAPQ >= 1 baseline placed (wrong) |
|---|---|---|---|---|---|
| simulation, rounds 23-25 | 12,585 (8,269) | 1,104 (0%) | 0 | 0 | 4,316 (0.14%) |
| HG002, donor truth | 27,199 (17,247) | 1,381 (1.09%) | 2 | 9 (8) | 8,958 (27.1%) |
| A119b fragments | 2,394 (1,108) | 319 (0.63%) | 4 | 5 (5) | 1,286 (13.8%) |

**Reading:**
- **O2 does not resolve exact ties (MAPQ 0).** Nothing in those reads separates the copies (0 / 2 / 4 placed), so they stay
  sets.
- **On near-ties O2 rarely overrules minimap2.** Its value is precision: on real reads, minimap2's MAPQ >= 1 near-tie calls
  are 14-27% wrong, and the subset O2 validates is about 1% wrong.
- **In one genome minimap2 is already right** (0.14% wrong), so O2 adds only certification and sets there.

The bars that would test this claim:
- held out on a donor never used here (the public testis library, ERR13885926);
- a minimum number of placements;
- O2's precision above the MAPQ >= 1 baseline's on the same reads;
- source reads counted, pseudoautosomal families apart.
