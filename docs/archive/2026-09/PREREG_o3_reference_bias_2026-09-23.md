# Pre-registration — O3: do copies absent from the reference hide as unmapped or multi-mapping reads?

**Written 2026-09-23 (§6zd), before any arm is run.** Advisor's claim: copy numbers match only when the sample
is the individual the reference came from; otherwise extra copies are unmapped or disguised as multi-mapping
reads. User's experience: SNP-level differences, not whole copies. Item 2 of `docs/PENDING_2026-09-23.md`.

## Arm A — simulation with truth (this document's decision arm)

**Sample genome:** human CHM13 chr20 plus, for each of K = 40 multi-exon RefSeq genes chosen at random among
those with ≥ 3 exons and a single annotated locus on chr20, one EXTRA COPY of every transcript of the gene,
mutated at divergence d ∈ {0.5%, 1%, 2%, 3%, 5%} (10 genes per level would confound gene with level, so
every gene gets a copy at EVERY level, in separate simulations): substitutions at random positions, no
indels (the question is identity, not structure). The extra copy is a transcript-level object: its reads
are the mutated transcript sequence, so mapping it to the unmodified chr20 is exactly the situation the
claim describes — a copy that exists in the sample and not in the reference.

**Reads:** `bench/sim_reads.py` HiFi-like errors (0.001, indel 0.0003), ends jittered ± 30 bp; 10 reads per
template transcript and 10 per extra-copy transcript. Mapped to the UNMODIFIED chr20 with the shipped
minimap2 settings (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`).

**Measured, per extra-copy read, per divergence level:** unmapped; mapped primary at the template locus
(absorbed) and its `de`; mapped primary elsewhere; MAPQ 0 / AS-tied (has a secondary within 0.98 of its AS);
and for the template reads the same, as the control. Then per gene: does the template locus show the
O3 signature — the extra-copy reads' `de` distribution separated from the template reads' (the divergence
rate O3's flag pass thresholds).

## The bar — committed now

| outcome at d ≤ 2% (the within-species range) | verdict |
|---|---|
| ≥ 80% of extra-copy reads absorbed at the template locus as primaries with `de` ≈ d, ≤ 10% unmapped, ≤ 10% tied | ⛔ **CLAIM REFUTED** — a missing copy is a divergence signal at the paralogue, which is O3's detector, not a hole in the data |
| 20-50% unmapped or tied | ⚠ **PARTLY** — the copy is visible but a fraction of its evidence is lost |
| ≥ 50% unmapped or tied | ⭐ **CLAIM SUPPORTED** — the pipeline must recover copies from unmapped/tied reads |

At d = 3-5% (cross-species-like) the same table is reported as information, not judged.

**Predicted, before looking:** ⛔ — minimap2 `splice:hq` aligns a 2%-diverged full-length transcript to its
template with MAPQ 60 and `de` ≈ 0.02 (the r1057 dumps show `de` up to 0.0135 on real MAPQ-60 reads); reads
go unmapped only when no 15-mer seeds survive, which needs divergence well above 5%; a tie needs a second
locus of equal score, which a single-copy template does not have. The extra copy therefore shows up as an
elevated-divergence pile at the template — exactly the O3 flag-pass statistic.

## Arms B and C (registered, not yet run)

- **B, same individual vs different individual, real data:** ⚠ corrected 09-23 — the fibroblast library
  (`fibroblasts/GCA_029281585.2_flnc_mm.bam`, the assembly's own animal KB3781) is ALREADY on the same GCF
  primary as the OR6737 testis library (`winloci_data/GGO_mm.bam`), so the same-vs-different-individual
  contrast needs no new genome-wide alignment (unmapped %, MAPQ-0/tied %, O3 candidates, copies per family;
  confounded by tissue). The DIPLOID TRUTH for the same animal exists on disk: `o3_hapcnv/{pat,mat}.chr.fa`
  (mGorGor1 v2.0 parental haplotypes) — reads from a copy present only on the haplotype the primary did not
  take (`pri_provenance.tsv`: 16 PAT / 9 MAT chromosomes) are genuinely reference-absent on the primary,
  and their pat/mat placement says where they belong. That is the real-data version of arm A.
- **C, cross-species:** A119b reads on the gorilla reference and vice versa (judge arm for any rule that
  arm A motivates).

---

# ARM A — OUTCOME (2026-09-23): ⛔ **claim refuted at every divergence tested; the extra copy is a divergence pile at its template**

chr20, 40 genes (106 transcripts), 10 reads per template transcript and 10 per extra-copy transcript per
level, shipped minimap2 against the unmodified chr20 (`/mnt/linuxdisk/tmp/gw22/o3/simA.py`; consolidated with simB/simC as `bench/o3_sim_copies.py`, mode `transcript`):

| divergence of the extra copy | extra-copy reads | unmapped | primary elsewhere | AS-tied / MAPQ 0 | **absorbed at the template locus** | median `de` at the locus (template reads: 0.0015) |
|---|---|---|---|---|---|---|
| 0.5% | 1,060 | 0 | 0 | 0 | **100%** | 0.0065 |
| 1% | 1,060 | 0 | 0 | 0 | **100%** | 0.0115 |
| 2% | 1,060 | 0 | 0 | 0 | **100%** | 0.0214 |
| 3% | 1,060 | 0 | 0 | 0 | **100%** | 0.0314 |
| 5% | 1,060 | 0 | 0 | 0 | **100%** | 0.0512 |

Not one read of a copy absent from the reference is lost or disguised, up to 5% divergence (beyond the
within-species range and into cross-species territory): every read maps as a MAPQ-60 primary to the
template, carrying its divergence in `de` almost exactly (0.0065 → 0.0512 for 0.5% → 5%, against 0.0015
for the template's own reads). That is the O3 flag-pass statistic — a per-locus excess of divergence — and
it is where a missing copy shows up. **⛔ by the bar (≥ 80% absorbed, ≤ 10% unmapped, ≤ 10% tied: measured
100 / 0 / 0).** Prediction confirmed.

Scope of what this settles: the mechanism claim for a copy with the SAME exon structure as its template (a
young duplicate). Two things it does not cover, registered for arms B/C: a copy whose template is itself
multi-copy (the tie is then between paralogues, not with "nowhere" — reads are still placed, at MAPQ 0 among
the family, which is O2's subject, not a loss), and a copy with structural differences (its reads would be
`j`-class chains at the template, also not unmapped). Unmapped reads need no seed to survive at all, which
is far beyond 5% divergence for full-length transcripts.

**⚠ Reconciliation with the 08-14 whole-genome excision panel (register rows 21/24, `project_o3_excision_wholegenome`):**
that experiment deleted one copy of 162 real two-copy gorilla families and realigned the matched fibroblast
reads: **64.2% of families ABSORBED** their orphans onto a paralogue (median concentration 0.967, depth
ghost 1.75×) but **33.3% ORPHANED** them (median 92.7% of the copy's reads became UNMAPPED). So the
advisor's "unmapped" mechanism is real when the nearest surviving relative is a distant paralogue (real
paralogues sit at median identity ~0.82, far beyond arm A's 5%); arm A's 100% absorption is the regime
of a young, structurally identical copy. The two together give the rule: **a missing copy's reads are
absorbed (with a `de` excess) below roughly 5-10% divergence to the nearest reference copy and unmapped
above it** — and the real BAMs carry almost no unmapped long reads (fibroblast 959 reads at median 69 bp;
testis 5,519 = 0.13%), so on this animal and on OR6737 the orphaned regime is nearly empty.
