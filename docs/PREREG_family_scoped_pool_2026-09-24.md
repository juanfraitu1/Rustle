# Pre-registration — family-wise assembly as a SCOPED secondary pool (§6zl), written before any arm is assembled

Advisor (via the user): the method should share information across the family / assemble family-wise. What is
already measured: admitting EVERY secondary into assembly (`RUSTLE_GTF_SECONDARY=1`, §6n2) makes members more
COMPLETE (NPIP complete chains 10 → 20/26) at a 9× transcript cost and, on gorilla NC_073244.2, doubles the share
of genes split into ≥ 2 loci (27% → 49%) and drops pair precision (1.000 → 0.903, row 1100's band scorer); the
tied pool (`AS ≥ 0.98 × genome-wide best`, r1060) finds more MEMBERS at precision 1.000 (≥ 90% identity pairs
18 → 36/108) but changes chains little. The untested middle: **seed with the tied pool, then pool ALL secondaries
only INSIDE an already-defined family, only for completing its members' isoforms.**

## Arms (gorilla NC_073244.2, `GGO_mm.bam`; the de novo chain of `docs/PREREG_locus_read_pool_2026-09-22.md`)

| arm | read pool for assembly |
|---|---|
| P | primaries |
| GOOD | primaries + secondaries with AS ≥ 0.98 × genome-wide best (existing `ggo44_GOOD0.98`) |
| **SCOPED** | GOOD's pool + every secondary whose placement overlaps a member locus of a GOOD family AND whose molecule's PRIMARY overlaps a member locus of the SAME family (families = `ggo44_GOOD0.98.fam.clusters.tsv`, size ≥ 2) |
| ALL | primaries + every secondary (existing `ggo44_ALL`, the ceiling and the cost) |

SCOPED is built as a filtered BAM (`ggo44_scoped.bam`) and assembled with `RUSTLE_GTF_SECONDARY=1` and no ratio
(every record of the filtered BAM is admitted); everything else is the chain's shipped polish. A secondary from
a molecule whose primary lies outside the family is an ECHO of a non-member and stays excluded — that is the
whole difference from ALL.

## Metrics — committed now

Universe: expressed referee genes (≥ 2 primary reads) that are members of a GOOD family (fixed set, so it cannot
move with the arm). gffcompare against `ref/NC_073244.2.ref.gtf`:
- **complete members** (primary) = universe genes with ≥ 1 transcript of class `=` (exact intron chain);
- partial = class `=`, `c` or `k`;
- cost = transcripts per member locus (SCOPED / GOOD), gffcompare transcript-level precision, loci per expressed
  referee gene (fragmentation, ≥ 2 loci share);
- pair precision vs the complete referee and ≥ 90% recall from `referee_band_score.py` (must not fall).

| outcome | verdict |
|---|---|
| complete members ≥ 1.2 × GOOD, transcripts per member locus ≤ 2 × GOOD, pair precision ≥ 0.98, ≥ 2-loci share ≤ GOOD + 5 points | ⭐ **family-wise completion works and is affordable: ship as the `families`-aware assembly option** |
| complete members ≥ 1.2 × GOOD but any cost bound fails | ⚠ **it completes members at ALL's price; report both, no default** |
| complete members < 1.2 × GOOD | ⛔ **within-family pooling does not complete members; the 5′ loss is the library's** |

**Predicted:** ⚠ — completeness rises (the §6n2 mechanism is exactly within-family echoes of full-length molecules)
but the transcript cost inside families is the same 9× per member locus, because the echoes ARE the completions.
Development substrate (human chr20, A119b) is run second with the same recipe if the gorilla arm is ⭐ or ⚠.

**Secondary outcome, added before SCOPED is assembled (the diagnosis below is of GOOD, not of SCOPED):** on GOOD the
referee pairs at 80–90% annotated-mRNA identity are lost 10% to a missing locus, **52% to an UNCLUSTERED locus (no
edge)** and 27% to a different cluster (317 pairs; ≥ 90%: 27 / 34 / 7 of 108). An unclustered locus with a partner
at 80–90% identity is the coverage clause failing on an incomplete representative — so if SCOPED completes members,
edge recall should rise. Pre-registered: **≥ 90% and 80–90% pair recall (`referee_band_score.py`) up by ≥ 0.05
absolute at pair precision ≥ 0.98 = a real family-definition gain from completeness**, reported alongside the
primary bar; below 0.05 = completeness does not reach the edges.

I will not change the pools, the universe, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-24) — `sec/ggo44_SCOPED.*` (39 s; filtered BAM = 159,348 primaries + 2,066 tied + 19,911 family-scoped secondaries over 28 GOOD families / 91 member loci)

| arm | complete members (=) of 73 | partial (=,c,k) | transcripts per member gene | gffcompare tx sens / prec | loci per expressed gene 0/1/2/3+ | ≥ 90% pair recall | 80–90% | pair precision |
|---|---|---|---|---|---|---|---|---|
| P | 55 | 57 | 5.58 | 26.5 / 36.8 | 113/365/107/28 | 18/108 | 20/317 | 1.000 |
| GOOD | 58 | 60 | 5.45 | 27.0 / 37.0 | 111/365/109/28 | 36/108 | 20/317 | 1.000 |
| **SCOPED** | **58** | 59 | 5.76 | 26.9 / 40.3 | 111/348/114/**40** | **28/108** | **15/317** | 1.000 |
| ALL | 58 | 60 | 27.01 | 27.3 / 19.8 | 79/271/130/133 | 38/108 | 30/317 | 0.903 |

## Verdict — ⛔ (complete members 58 = GOOD, bar was ≥ 70; secondary outcome: band recall FALLS, ≥ 90% 36 → 28)

The prediction (⚠, completeness up at ALL's transcript price) was wrong on both halves. Family members on this
contig are already 79% complete under GOOD (58/73 exact intron chains) and even ALL, the ceiling, adds none: the
§6n2 effect (NPIP 10 → 20/26 complete chains) is a property of that human family's 5′-truncated library, not a
general lever. What the within-family pool does instead is FRAGMENT: genes split into ≥ 3 loci rise 28 → 40 (the
echoes of a member's molecules at its paralogues become extra loci inside the family), `collapsed` 12 → 18, and the
representative each locus contributes changes, so eight ≥ 90% pairs and five 80–90% pairs that GOOD joined are no
longer joined. Transcript precision rises 37.0 → 40.3 only because the extra loci carry annotated echoes.
**Family-wise pooling of secondaries into assembly is closed on this substrate; the completion it was meant to buy
is already there.** The 52% of 80–90% pairs lost to an unclustered locus (r1101) are therefore not incomplete
representatives waiting for more reads — the next diagnosis is what those representatives' alignments actually
fail on (identity, coverage, or no record), which is the chr16 miss-diagnosis (r1096) repeated on gorilla.

## Post hoc — what the unclustered loci fail on (GOOD arm, DNA-level hull all-vs-all, not a rule change)

Of the 80–90% referee pairs whose two genes both have a locus but one joins no cluster (177): **140 (79%) have an
alignment record at ≥ 0.80 identity that covers < 50% of the shorter locus hull**, 16 have no record, 17 a record
below 0.80, 4 a full-coverage record dropped by the exon clause / MCL (e.g. LOC101133546–ZNF813: identity 1.0 over
60.6 kb, coverage 0.70). ≥ 90%: 21 / 9 / 2 / 6 of 38. The same shape as chr16's direct-edge misses (r1096):
partial alignments between loci of unequal extent, not seeding and not completeness. The coverage floor itself
was already searched (r1053/r1054: every denominator direction; r1097: union coverage at 0.30 over-merges the
held-out catalog; the shipped `--min-cov-shorter 0.70` moves the other way) — this is recorded as the diagnosis,
not as a proposal.
