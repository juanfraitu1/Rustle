# Where rustle+VG beats StringTie via multi-mapping reads — genome-wide portfolio

**2026-06-05.** Goal: find paralog loci where rustle's VG/multi-mapping machinery recovers
real transcripts that StringTie (which drops secondary alignments) structurally misses.

## Method
1. **Candidates**: regenerated the `expressed_real_copy` set from the paralog secondary-dependence
   scan (`bench/paralog_secondary_scan/scan_out/`) — 93 annotated multi-exon genes expressed
   predominantly via secondary alignments where the read sequence favors *this* copy.
2. **Per-locus eval** (`winloci_eval.sh`): extract candidate (+dominant sibling) region from the
   genome BAM; run StringTie `-L`, rustle baseline `-L` (primary-only), and rustle VG
   (`--vg --vg-snp RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1`); score each vs RefSeq
   (`GGO_genomic.gff`, class `=`).
3. **Read-backing verify** (`winloci_verify.sh`): strict edit-distance decisiveness — how many reads'
   sequence STRICTLY prefers this copy over siblings (the ZFY-real vs DAZ3-phantom discriminant).
4. **Four-filter portfolio** (`winloci_aggregate.py`): a genuine win transcript must be
   (a) VG-specific (in VG output, absent from rustle primary-only baseline),
   (b) annotated at the CANDIDATE (secondary-dependent) locus,
   (c) NOT matched by whole-genome StringTie (not a region artifact),
   (d) read-backed (verify verdict = real).

## Authoritative sweep (deterministic, direct — NOT agent-wrapped) — 93 loci
| class | count |
|---|---|
| win_vs_st (raw) | 20 |
| tie | 56 |
| regress (VG loses real tx ST keeps) | 16 |
| win_vs_base_only | 1 |

Raw transcript accounting vs StringTie: **+36 gained** (across 20 win loci) / **−29 lost**
(across 16 regress loci) = **net +7**, but with heavy churn and many gains that fail the
rigorous filters.

## FINAL PORTFOLIO — 4 loci, 7 VG-specific read-backed transcripts genome-ST misses
| gene | family / sibling | owner_frac | strict-owner reads | transcripts |
|---|---|---|---|---|
| **RABL2A** | RABL2 / RABL2B | 0.560 | 28 | XM_055379901.1, XM_055379919.2, XM_055379922.2 |
| **LOC129523511** | RABL2 / RABL2B | 0.409 | 27 | XM_055384089.2, XM_055384090.2 |
| LOC115933307 | / LOC101143875 | 0.133 | 11 | XM_055351714.2 |
| ZC3H7A | / LOC115932749 | 0.120 | 41 | XM_063699488.1 |

**Flagship = the RABL2 family.** RABL2A: VG matches 5 RefSeq tx vs StringTie's 4 (+1 net), recovering
3 transcripts ST misses; 28/50 chain reads strictly anchor to this copy. rustle baseline (primary-only)
also misses these 3 — so the **secondary/multi-mapping reads are doing the work** (the whole thesis).

## Honest caveats (the result is narrow and costly)
- **VG is NOT a free win.** It regresses at 16 loci (loses ~29 real tx, e.g. the RBMY array copies
  each lose 5) — more loci than it cleanly wins. The genuine value is *targeted recovery of
  identifiable paralog copies*, not a general sensitivity gain.
- **Window-fragile.** The win is sensitive to the exact extraction window: on the RBMY array, shifting
  the right edge ~20 kb flips a +1 win into a −5 regression. RBMY is therefore a poor flagship and
  drops out of the strict portfolio entirely. This is the documented productionization gap.
- **Region regime, not whole-genome.** These wins are measured on isolated candidate regions (the
  favorable regime where tandem/family resolution fires). Whole-genome VG behavior (bundles
  mis-grouped into unrelated families) is the open productionization question.
- **Methodology note:** an earlier LLM-agent-wrapped screen reported 30 wins; deterministic direct
  re-runs showed it over-reported (agents deviated from the exact command). rustle itself is
  deterministic. Lesson: run deterministic bash directly, not wrapped in agents.

Artifacts: `PORTFOLIO.json`, `sweep_summary.tsv`, scripts in `bench/multi_copy_eval/winloci_*.sh|py`.
Genome assets staged on ext4 at `/home/juanfra/winloci_scratch/` (9p `/mnt/c` is too slow for the FASTA load).
