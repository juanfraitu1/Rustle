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
