# V4c — identity-gradient frontier demo

A synthetic, ground-truth demo of Rustle's copy-resolution mechanism across the *full*
pairwise-divergence range, from clearly-separable copies down to byte-identical ones. It
answers the question an advisor asks first: **"does copy-number detection survive when
sequence identity gets high enough that per-read assignment can't?"**

Everything here is real output from the actual `copy_assign` binary and `minimap2` — no
numbers are hand-typed or interpolated. Every run is seeded (`random.Random` only, no
wall-clock), so re-running reproduces the same story (small individual-read jitter aside,
see **Caveats** below).

## What it shows

**Component A — the frontier curve.** A synthetic 4-copy tandem gene family (each copy:
`exon1(700bp) + intron(200bp) + exon2(700bp)`, same relative layout, own genomic
coordinates — this is what lets minimap2's splice-aware alignment tell the copies apart at
all; see `sim_tandem.py`'s module docstring for why an intron-less design fails) is built
at six target pairwise identities (~90/94/96/98/99.5/100%). At each point: simulate HiFi-like
reads (40/copy, 0.3% per-base error), align with `minimap2 -ax splice:hq -N 50`, run
`copy_assign --homology-primary --min-copies 2 --dump-psv` on the *intact* reference (no
deletion — ordinary assignment), and record what came out of `families.tsv` (copy number)
and `assignments.tsv` (assigned vs. tied).

**The result** (`bench/mechanism/verification_results.json["v4c"]["ladder"]`):

| achieved identity | copies detected (planted=4) | fraction assigned | fraction tied (abstained) |
|---:|:---:|---:|---:|
| 90.1% | 4 ✓ | 100.0% | 0.0% |
| 94.0% | 4 ✓ | 100.0% | 0.0% |
| 96.0% | 4 ✓ | 99.4% | 0.6% |
| 98.0% | 4 ✓ | 26.0% | 74.0% |
| 99.6% | 4 ✓ | 25.0% | 75.0% |
| 100.0% | 4 ✓ | 0.0% | 100.0% |

Copy **number** is recovered correctly (4/4) at *every* point on the ladder, including at
100% identity where the copies are byte-identical — this comes from read-conflict topology
(which loci a read's multiple equally-good alignments span), not from sequence content, so
it survives even where sequence content carries zero information. Per-read **assignment**
degrades sharply once pairwise identity crosses roughly 96-98%: fewer and fewer reads span
a PSV column or junction that actually differs between copies, so the certificate correctly
abstains (`status=tied`) rather than guess. At 100% identity every read abstains — the
honest floor, not a bug. This two-line divergence (flat green, rising red) *is* the
identifiability frontier, and it is also the built-in rebuttal to "your method just
collapses everything together": copy number stays right even when assignment gives up.

See `frontier.svg` for the plotted version (x = achieved pairwise identity, y = fraction;
green = copy-number-correct, red = fraction of reads abstained/Tied).

**Component B — the recovery anchor.** Proves the *other* half of the story: that a copy
deleted from the reference is genuinely *recovered* (not just detected as "something's
there"), at a divergence where recovery is possible. This reuses `sim_tandem.py`'s
already-proven V4a recipe verbatim (a fresh subprocess run, not cached numbers): a 4-member
family at ~2% per-copy divergence (~96.5% mean pairwise identity), of which 2 members stay
in the reference (`A`, `B_host`, forming a real 2-copy family) and 2 are deleted
(`hidden1`, `hidden2`) and pool their reads onto `B_host`'s locus. Run with
`copy_assign --min-copies 2 --absent-copies --linearize --dump-psv --phase`. Both hidden
copies are **admitted** (not merely flagged as "DNA needed") and their PSV-column alleles
match their true deleted sequences at identity **1.0** (`set_match`, permutation-invariant).

```json
"recovery_anchor": {
  "identity": 0.9652,
  "recovered": true,
  "min_identity": 1.0,
  "n_psv_cols": 96
}
```

### An honest, disclosed detour: why the anchor isn't literally "the most divergent ladder point"

The brief's first instinct (put the recovery anchor at the ladder's *lowest* identity,
~90%, since that "sounds like" the hardest/most-interesting case) turned out to be
mechanically backwards, and this demo keeps that finding rather than hiding it
(`verification_results.json["v4c"]["recovery_anchor"]["probes_tried_first"]`):

- **Probe 1** (~90% identity, the ladder's most divergent point): `copy_assign` logged
  `conflict-graph (de-tie): 0 edges -> 0 families`. The admission/recovery certificate
  attaches a reference-absent candidate to an *already-detected* family, and that family
  is formed here by the plain read-conflict graph, which needs residual MAPQ-0 *ambiguity*
  between loci. At ~90% identity every read resolves cleanly to its true locus — there is
  no ambiguity left to build a family from, so there is nothing to attach a recovered copy
  to. High divergence makes recovery unnecessary (assignment already works) and, in this
  design, impossible to trigger.
- **Probe 2** (retargeted to ~98%, this demo's own analytical host+2-hidden partition
  search over the ground-truth identity matrix — see `_pick_host_hidden_partition` in
  `sim_ladder.py`): logged `3 skeletons -> 3 transcripts -> 1 reps` — the searched
  partition's two "kept" loci collapsed into a single assembled representative instead of
  staying distinct, so again 0 conflict edges. A real gap in that search heuristic (it
  optimizes how tightly the hidden pair clusters onto the host, but does not independently
  guarantee the host and the other anchor stay far enough apart to remain two distinct
  loci), disclosed rather than patched around under time pressure.
- **What actually ran**: V4a's own hand-vetted recipe (`sim_tandem.py::_divergent_arm`,
  `build_tandem(seed_seq, 4, 0.02, seed_rng=11)`), which already empirically produces a
  clean 2-locus split with the hidden pair pooling 38/40 and 39/40 onto the host — reused
  here as a **fresh** subprocess run (not a copy of old numbers), landing at 96.5% mean
  pairwise identity — inside this ladder's own 96-98% band, i.e. exactly where the frontier
  curve above shows assignment starting to abstain. That is the mechanically right place
  for a recovery demo to live: enough residual ambiguity to form a family and admit a
  candidate, but real divergence to distinguish what's recovered from what it's
  confused with.

## Files

- `run_demo.sh` — regenerates everything (ladder + recovery anchor + JSON + SVG + example
  inputs) from scratch. `./run_demo.sh` for the full 6-point ladder, `./run_demo.sh --fast`
  for a 4-point subset (94/98/99.5/100%) if the full sweep is too slow in your environment
  (both are already fast in practice — the full 6-point ladder + recovery anchor runs in
  well under a minute; `--fast` exists mainly for constrained CI-like environments).
- `frontier.svg` — the two-series plot described above, hand-emitted plain SVG (no
  plotting library).
- `ref.fa` / `reads.fq` — the example inputs for the ladder's most-divergent point
  (~90% achieved identity, 4 copies × 40 reads = 160 reads, ~7KB reference), committed so
  the demo is inspectable without re-running anything. Regenerated fresh by `run_demo.sh`
  on every run (same seeds ⇒ same bytes).
- `../sim_ladder.py` — the driver (imports `build_tandem`/`set_match`/etc. from
  `../sim_tandem.py`, the V4a module; does not reinvent that machinery).

## How to reproduce a single ladder point by hand

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
/home/juanfra/miniforge3/bin/minimap2 -ax splice:hq -N 50 \
    bench/mechanism/demo/ref.fa bench/mechanism/demo/reads.fq \
  | /home/juanfra/miniforge3/bin/samtools sort -o /tmp/x.sorted.bam -
/home/juanfra/miniforge3/bin/samtools index /tmp/x.sorted.bam
target/release/copy_assign \
  --bam /tmp/x.sorted.bam --fasta bench/mechanism/demo/ref.fa \
  --region ladder0:1-7000 --out /tmp/x_out \
  --homology-primary --min-copies 2 --dump-psv
column -t /tmp/x_out.families.tsv
cut -f4 /tmp/x_out.assignments.tsv | sort | uniq -c   # status breakdown
```

(Substitute `/tmp/x.sorted.bam` etc. for a path under `~/winloci_scratch/` if you'd rather
follow this repo's WSL2 crash-avoidance convention of never writing scratch output to
`/tmp`; the demo driver itself always does.)

## Caveats

1. **`n_reads` in the ladder table counts assignment *rows* (one per read × alignment
   record), not distinct read molecules.** At high identity, `minimap2 -N 50` reports up
   to 4 near-tied alignment records per physical read (one per equally-good locus), and
   `copy_assign` emits one `assignments.tsv` row per record — so `n_reads` climbs from 160
   (all-primary, low identity) to 640 (4 records × 160 reads, 100% identity) even though
   the same 160 simulated reads (`n_distinct_reads`, always 160) are the only physical
   input. This is itself part of the story (more equally-good placements = more
   ambiguity), not an artifact to hide, but it means `frac_assigned`/`frac_tied` are
   fractions of *records*, and that is disclosed here rather than left implicit.
2. **Minor read-level nondeterminism.** Re-running the full ladder twice gave qualitatively
   identical curves but occasionally flipped a single read's status at a tie boundary
   (e.g. 99.4% vs. 100.0% assigned at the 96%-identity point across two runs). This traces
   to `minimap2`'s own multi-threaded chaining being non-deterministic on exactly-tied
   alignment scores, not to anything in `copy_assign` or this driver (both of which are
   pure/deterministic given a fixed BAM). The shape of the frontier is robust to this; the
   exact fourth decimal place of any one point is not guaranteed bit-identical run to run.
3. **The recovery anchor is not literally "the most divergent ladder point"** — see the
   dedicated section above. This is a real, source-confirmed constraint
   (`min_clusters=3`, hardcoded in `absent_copy.rs`, requiring residual read-conflict
   ambiguity to even form the family a candidate attaches to) inherited from V4a, not new
   to this demo, but worth restating: recovery-via-certificate lives in the same
   high-identity band where assignment starts to abstain, not in the easily-resolved
   low-identity band.
4. Ladder identities are reported as **achieved**, not target (bisection search gets close
   but not exact — e.g. target 0.90 achieved 0.9011). Never assume; every identity number
   here is measured by counting differing bases between the planted ground-truth copy
   sequences (`_mean_pairwise_identity` in `sim_ladder.py`), which is exact (no alignment
   noise) since the true sequences are known by construction.
