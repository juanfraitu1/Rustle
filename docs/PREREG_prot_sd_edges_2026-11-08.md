# PREREG row 1590: protein edges between SD-linked loci on top of SD damping, held out on HG002

Written 2026-11-08, before any binary run or locus-level protein table exists. The setting is chosen on A119b (dev) by
the rule below and appended here before the held-out jobs start. Judged against row 1580's damped curve, like row 1581.

## Question

On top of `--sd-damp 0.25` (row 1580, PASS on HG002), does an edge between two loci that (a) one SD pair links and
(b) whose proteins align buy recall more cheaply than lowering MCL inflation does? Row 1581 added an edge between EVERY
SD-linked pair without a homology edge and failed on HG002 (-0.007); the protein condition is the new ingredient.

## Provenance (row 1589, step 0, post hoc)

Protein links alone: 206 of 565 missed truth pairs reached, but 5,907 new pairs at precision 0.035 (old paralogs),
STOP as registered. Restricted post hoc to gene pairs one CHM13 SEDEF pair holds: 103 reached, 169 new pairs, precision
0.609 (direct pairs, CAT gene exons, no MCL). This row tests the same idea on de novo loci, inside MCL, held out.

## Rule (`bench/denovo_eval/protedge.py`, then the existing `mcl_families --extra-edges`)

1. Locus proteins: every transcript of the run's `families.gtf` goes to the locus of `PREFIX.fam.loci.gff3` (same
   strand) sharing the most exon bases with it; its longest ATG..stop ORF >= 80 aa on its strand (an ORF running off
   the 3' end counts) is kept; <= 3 distinct ORFs per locus, longest first.
2. DIAMOND 2.2.8 `blastp --more-sensitive -e 1e-10 -k 0`, all vs all. Two loci are protein-linked when an ORF pair has
   e <= 1e-10, alignment length >= 50% of both proteins and identity >= 50%.
3. SD gate: kept only if one SEDEF pair links the two loci by `SdPairs::holds` (one copy holds >= half of a locus's exon
   bases AND >= half of its span, the other copy the other locus), computed on the gff3 exon blocks. Check before any
   arm is scored: this port's SD-linked pairs among the gated-PAF edges must be within 1% of the binary's
   `sd_damp_edges_linked` on the same run, else stop and fix the port.
4. Each kept pair -> one line `locusA locusB W` for `--extra-edges` (an existing edge is raised to max(edge, W)).
   Extra edges are merged before `--sd-damp`, which leaves SD-linked edges undamped. Gate: the binary must report 0 new
   nodes (every name matched).

## Setting (fixed on A119b before any held-out run)

Candidates: W in {0.25, 0.5, 1.0} x I in {2.8, 3.2, 3.6}, all with `--sd-damp 0.25` and the CHM13 SEDEF BED. Baseline:
the binary damped-only curve (r1578span / r1581dev: I = 1.6 .883/.736, 2.0 .909/.716, 2.4 .914/.722, 2.8 .925/.718,
3.2 .943/.716, 3.6 .943/.688), `pairscore.py` on PREFIX.fam.copies.tsv. Choose the candidate with the largest
precision gain over that curve at its recall (linear interpolation between the two bracketing points, points sorted by
recall) among candidates with recall >= .718 (the damped default). Stop, no held-out run, if no candidate's recall lies
inside the curve's recall range or the best gain is < +0.04. Row 1581 went to held-out with its dev gain equal to the
bar and failed; this row asks for a margin.

## Held-out

HG002 (decisive): `r1580/hg002` nodes and the HG002 SEDEF build (`sedef_build/hg002pat/final.bed`, 4 chr9-satellite
candidates dropped), the chosen W at I in {1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0}; baseline = the HG002 damped
curve already run for row 1580 (`r1580/hg002/sd_I*`). Gorilla OR6737 descriptive only (88 truth pairs).

## Bar (HG002)

PASS when both hold:

1. Precision of the chosen arm >= the HG002 damped curve's precision at its recall (bracketing interpolation) + 0.03.
2. Recall of the chosen arm >= recall of the HG002 damped default (`sd_I2.8`, .321).

If its recall lies outside the damped curve's recall range: FAIL.

## What follows

- PASS: protein-SD edges documented as an opt-in step (`protedge.py` + `--extra-edges`); a binary port and defaults
  are the user's call.
- FAIL: register; the script stays in bench as a recorded negative.

## Outcome (2026-11-08, A119b dev, jobs 4453782 / 4453783): STOP, no held-out run

Port check (on the binary's own graph, `--dump-graph` of d25 I2.8): this port finds **115,510** SD-linked edges of
525,890, the binary's `sd_damp_edges_linked` is 115,510 (exact). 531,400 protein-linked locus pairs, 50,936 of them
SD-linked (48,453 already a PAF edge): 3,618 edges added per arm, 0 / 869 / 44,705 raised at W 0.25 / 0.5 / 1.0.
**Node gate FAILED as written:** 154 new nodes, not 0. All 154 are loci of this run's `fam.loci.gff3` that had no
homology edge (graph nodes exist only with an edge), so no name failed to match (the gate's purpose), but those loci join
families through protein edges alone. The port check ran in a job parallel to the arms, not before them. Neither changes
the verdict below, which is a stop either way.

| arm | precision | recall | damped curve at recall | gain |
|---|---|---|---|---|
| W0.25 I2.8 | .918 | .729 | .899 | +.020 |
| **W0.25 I3.2** | **.936** | **.728** | **.901** | **+.035** |
| W0.5 I2.8 | .897 | .735 | .885 | +.012 |
| W0.5 I3.2 | .928 | .725 | .907 | +.021 |
| W1.0 I3.6 | .913 | .732 | .892 | +.021 |
| W0.25 I3.6, W0.5 I3.6 | .935, .929 | .709, .716 | | recall < .718, not eligible |
| W1.0 I2.8, W1.0 I3.2 | .891, .908 | .737, .743 | | recall above the curve's range (max .736) |

Best eligible gain +0.035 < the +0.04 dev margin: **STOP as registered**, no HG002 run. Post hoc, no claim: W1.0 I3.2
reaches recall .743 at .908, beyond every damped point (.736 at .883); the curve cannot price it.
