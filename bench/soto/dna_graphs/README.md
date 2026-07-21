# Soto DNA variation graphs -- all families

76/83 Soto families rendered as a DNA variation graph (base-level abpoa MSA of the members' genomic sequences -> GFA); 7 skipped (logged below and in `INDEX.tsv`, never silently dropped).

## What the colours mean

- **green** (`#1e8e3e`) -- a copy (member) that RNA (IsoSeq A119b) actually RECOVERED (`soto_member_detection.tsv`, `detected=Y`).
- **red** (`#d93025`) -- a copy that RNA did NOT recover (DNA-only / K=0-floor); its per-member cause (from `soto_floor_decomposition.tsv`) is in each family's `<stem>.legend.tsv`.
- **grey** (`#9aa0a6`) -- a shared/conserved backbone **node** in the graph that BOTH a green and a red copy pass through (mixed). Grey is a **node**-level colour only -- every member (copy) itself is always definitively green or red per RNA detection; a node is grey where green and red copies' sequences coincide.

The DNA variation graph is the identifiability CEILING: every Soto DNA-catalog copy is a path in the graph regardless of whether RNA recovered it. The graph REPRESENTS what is given (all copies, from Soto's independent DNA-read-depth catalog); it does not 'detect' families.

## How to open a GFA

Open any `bench/soto/dna_graphs/<family_id>_<gene>.gfa` in [Bandage](https://rrwick.github.io/Bandage/) (`File > Load graph`), then `File > Load CSV` the matching `.colours.csv` and colour nodes by the `Colour` column (Bandage: Graph drawing -> Node colours -> `Custom colours based on...` -> the loaded CSV column). The `.legend.tsv` alongside each GFA lists, per member: detected Y/N, recovered_by, cause (for red members), colour.

Look up any family fast in `INDEX.tsv` (`family_id`, `gene`, `n_members`, `n_green_rna_recovered`, `n_red_dna_only`, `n_grey`, `gfa_path`, `status`, `skip_reason`).

## Summary counts

- **76/83** families rendered, **7/83** skipped.
- Across rendered families: **355** members total -- **273 green** (RNA-recovered) / **82 red** (DNA-only) / **18465 grey graph-nodes** (shared/mixed backbone, node-level not member-level).
- These counts are over RENDERED families only -- members of a SKIPPED family (e.g. the 7 single-member families that can't form a comparative graph) are not pictured anywhere and are excluded from the green/red totals above; see `INDEX.tsv` `n_members` for every family including skipped ones, and `bench/soto/soto_member_detection.tsv` for the ground-truth green/red split across ALL 362 Soto members genome-wide (the 76.2% RNA-recovery figure).

## Skipped families

| family_id | gene | n_members | reason |
|---|---|---|---|
| ID_141 | AC119751 | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_142 | CR381670 | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_222 | AC243829 | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_267 | AL590399 | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_302 | BOLA2B | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_334 | DEFB104B | 1 | single-member family (n_present=1) -- no comparative graph to draw |
| ID_348 | DUX4L50 | 1 | single-member family (n_present=1) -- no comparative graph to draw |

## Regenerating

GFAs are regenerable (deterministic given the same inputs): `/home/juanfra/miniforge3/bin/python3 bench/soto/soto_dna_vg.py` (needs `pyabpoa`; run from the repo root). Each family renders in its own child process so one abpoa OOM/crash cannot take down the run.


## Reads-level graph (DNA ceiling vs what the reads resolve)

For **AMYLASE (ID_131, chr1:103.4–103.8 Mb, a tight 6-copy tandem array)** there is a matching
READS-level copy graph so you can show the DNA ceiling beside what IsoSeq actually resolves:
- `reads_ID_131_AMYLASE.copygraph.gfa` — the 13 copy-paths + PSV bubbles the reads support
  (the full per-read-walk graph was 756 MB / 33,744 read-walks; this compact version keeps the
  copy paths + PSV segments only, ~1.4 MB, Bandage-loadable) + its `.colours.csv`/`.legend.tsv`.
- `reads_ID_131_AMYLASE.exon.gfa` — the copy-specific exon graph (29 KB).
- Real run: `copy_assign --phase --dump-psv` on A119b IsoSeq vs CHM13, 1 family / 12 copies /
  12 haplotypes phased.

The other two flagships (SRGAP2 ID_462, PMS2P ID_8) are **dispersed** across their chromosome
(84 Mb / 33 Mb spans, copies at distinct far-apart loci — not co-located arrays), so a single
`copy_assign --phase` region is not meaningful for them; their DNA graphs above already represent
all copies. AMYLASE is the clean co-located example for the DNA-ceiling-vs-reads comparison.
