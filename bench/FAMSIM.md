# famsim — controlled gene-family simulations that prove the condition they test

**2026-10-02.** Design: `docs/FAMSIM_DESIGN_2026-10-02.md`. Code: `bench/famsim/` (Python 3, stdlib + pysam + numpy/scipy;
minimap2 and samtools on PATH; the Rust binaries from `--bin` / `RUSTLE_BIN`, default
`/mnt/linuxdisk/home/juanfraitu/rustle_target/release`). Tests: `python3 bench/famsim/test_famsim.py` (16, < 1 s).

The advisor's objection is that the results pick easy cases. His demonstration: plant a gene A and an identical copy A′
in an artificial chromosome, show what minimap2 does, then mutate A′ — SNVs, then more or fewer exons, inversions — and
show the pipeline still works, for any family and any ape. famsim does exactly that, and its `verify.tsv` proves from the
output files that each condition was really simulated, so minimap2's own limits are its problem, not a doubt about the
simulation.

## Quick start

```bash
python3 bench/famsim spec > my.json                      # a template spec to edit (NPIPA1 from CHM13 RefSeq, chr20 background)
python3 bench/famsim spec exon_loss > rung.json          # one rung of the ladder as a spec
python3 bench/famsim all my.json --out /mnt/linuxdisk/tmp/fs/my --threads 2      # make + verify + align + run + score
python3 bench/famsim ladder --template-spec bench/famsim/examples/ladder_human_SNRPB.json --out /mnt/linuxdisk/tmp/fs/human
```

`ladder` resolves the template gene and the decoys once (`template.json`, `decoys.json`), runs every rung in its own
directory and writes `ladder.tsv` (one row per rung: verify status + the headline numbers). `--only r1,r2` reruns rungs
in place. Each rung takes 1-3 s on a 300 kb contig; a 22-rung ladder ~30-60 s.

Individual steps: `make SPEC --out DIR`, `verify DIR`, `align DIR`, `run DIR [--stages assemble,denovo,guided,assign,flag]
[--copy-table families|catalog]`, `score DIR`.

## Scenario spec

```json
{"name": "exon_loss", "seed": 7,
 "background": {"source": "fasta", "path": ".../chr20.fa", "region": "chr20:20000000-20300000"},
 "template":   {"source": "annotation", "genome": ".../chm13v2.0.fa", "annotation": ".../x.gff.gz", "gene": "NPIPA1", "intron_cap": 3000},
 "copies": [{"id": "A"},
            {"id": "A2", "ops": [{"op": "snp", "rate": 0.02}, {"op": "exon_delete", "exon": "middle"}]},
            {"id": "A3", "contig": "sim2", "strand": "-", "in_reference": false, "expression": 0, "ops": []}],
 "decoys": {"n": 3, "min_exons": 3, "max_span": 30000, "chrom": "chr20"},
 "layout": {"spacing": 20000, "start": 20000},
 "reads": {"per_copy": 50, "err": 0.001, "indel": 0.0003, "jitter": 30, "trunc5_frac": 0.0, "trunc5_max": 0.3}}
```

- **background** `fasta` (a real slice; contig `sim2` takes the next slice of the same length) or `random` (`length`, `gc`).
  The builder refuses a slice that contains the real locus of the template or a decoy.
- **template** `annotation` (GFF3 RefSeq/CAT/Liftoff or GTF, plain or `.gz`; `gene` name, `transcript` id, or random with
  `chrom` + `min_exons`/`max_exons`/`max_span`; `intron_cap` shortens long introns keeping 20 bp ends), `synthetic`
  (`exons`, `introns` length lists, GT…AG introns), or `file` (a saved `template.json`). The canonical transcript = most
  exons, tie longest. The template is the span first exon → last exon in transcript orientation.
- **copies** `id`, `pos` (0-based; default: laid out in order at `layout.spacing`), `contig` (default `sim`), `strand`,
  `in_reference` (false = O3 case: the copy sits on its own contig that `genome.ref.fa` omits), `expression` (reads; 0 = the
  semi-guided case), `isoforms` (`[{"skip": [3], "weight": 0.3}]` adds exon-skipping isoforms), `ops`.
- **decoys** `n` random unrelated genes (same annotation/genome/chrom as the template unless given), planted alternately on
  both strands and expressed like copies — the precision control.
- **reads** per copy; HiFi substitutions/indels; MANDATORY end jitter (identical reads collapse under dedup, §6n0);
  optional 5′ truncation (`trunc5_frac` of reads lose up to `trunc5_max` of their length).

### Operators (applied in order; exons numbered 1..n in transcript order, or `"middle"`/`"first"`/`"last"`/a label)

| op | arguments | DNA | RNA |
|---|---|---|---|
| `snp` | `rate`, `region` all\|exons\|introns | round(rate·L) substitutions, GT/AG protected | PSVs |
| `indel` | `rate`, `max_len` | 1..max_len bp indels inside one exon/intron | PSV indels |
| `exon_delete` | `exon` | exon removed (terminal: + its intron); flanking introns merge GT…AG | one exon fewer |
| `splice_kill` | `exon` (internal) | donor GT→CT | exon skipped, DNA kept |
| `exon_insert` | `after`, `length` \| `source_exon`, `offset` | AG+exon+GT inserted into that intron | one exon more (`ins<k>` / `dup<label>`) |
| `exon_shuffle` | `a`, `b` | the two exon sequences swap places | exon order changes |
| `invert` | `exon` \| `intron` \| `span [a,b]` \| `whole` | reverse complement in place (intron: interior only, 6 bp ends kept) | exon/span: inverted exons leave the chain (`keep_in_rna` keeps their rc); intron: silent; whole: strand flip |
| `truncate` | `side` 5\|3, `exons` \| `bp` | that end removed (a cut inside an intron moves to the next exon) | partial copy |
| `convert` | `from` (an earlier copy), `exon` \| `span` | the donor's CURRENT sequence copied in | mosaic |
| `intron_resize` | `intron`, `length` | middle trimmed/filled, 20 bp ends kept | silent |

## Products of `make`

`genome.truth.fa` (every contig) · `genome.ref.fa` (+.fai; reference-absent contigs omitted) · `truth.gtf` (one
transcript per isoform, `exon_label`) · `truth.gff3` (gene/mRNA/exon with `Name=` and `gene=`, what `mcl_families --gff`
reads; the guided node set) · `truth.families.tsv` (`Gene Name / Family ID / Contig`, `family_score`'s format) ·
`copies.tsv` · `copies.fa` (each planted sequence, transcript orientation) · `reads.fq` (`copy|isoform|i`) ·
`reads.truth.tsv` (per read: chain interval, the junctions it contains) · `manifest.json` · `template.json` · `decoys.json`.

## `verify.tsv` — the proof

`claim  copy  expected  observed  PASS|FAIL|INFO  note`, measured on the products only:

| claim | measured how |
|---|---|
| `planted_in_genome` | span and every truth exon read back from `genome.truth.fa` at `truth.gtf` coordinates equal the copy's sequence |
| `junctions_canonical` | GT…AG at every truth intron on its strand (a template's own non-canonical motif is allowed) |
| `rna_exon_count`, `length_bookkeeping` | iso0 exon count / copy length = template + the op deltas |
| `pair_identity_genomic`, `pair_identity_chain` | `minimap2 -cx asm20 --eqx -X -N 50 -p 0.1 --secondary=yes` all-vs-all of `copies.fa` and of the spliced chains: nmatch/alen, gap-compressed 1−de, coverage; snp-only pairs must land within 35 % of the rate (+0.003); identical pairs at 1.0000; structural pairs INFO |
| `decoy_unrelated`, `background_clean` | no hit ≥ 300 bp at ≥ 0.80 identity between a decoy and A, or between a copy and the genome outside the planted intervals |
| `deleted_exon_absent`, `deleted_exon_not_in_rna`, `splice_killed`, `killed_exon_in_dna`, `inserted_exon`, `exons_shuffled`, `segment_inverted`, `inverted_exon_not_in_rna`, `truncated`, `exon_converted`, `whole_inverted` | 11-mer presence / orientation of the template exon in the copy, flank dinucleotides, label order in iso0, reverse complement read back from the genome |
| `absent_from_reference` | the copy's best hit matches fewer bases in `genome.ref.fa` than in `genome.truth.fa` (an identical absent copy is undetectable and FAILS honestly) |
| `read_count`, `read_junctions_in_truth` | reads per copy = expression; each read's junctions ⊆ its isoform's |

`all`/`ladder` stop a scenario at the first FAIL (`--keep-going` overrides).

## `run` — the shipped recipes

| stage | what runs | product |
|---|---|---|
| `assemble` | `tools/rustle_pipeline.sh assemble` (as_table + `copy_assign --assemble-only --genome-wide`, shipped polish, f1v2) | `run.gtf`, `run.families.gtf` |
| `denovo` | `tools/rustle_pipeline.sh families` (`mcl_families --from-gtf --emit-units`) | `run.fam.clusters.tsv`, `run.fam.copies.tsv/.fa` |
| `guided` | gene regions of `truth.gff3` → `samtools faidx` → `minimap2 -x asm20 -c --eqx -P` → `mcl_families --paf --gff` (the recipe of `figures/_o1_recovery.py::guided_families`) | `guided.clusters.tsv` |
| `assign` | `copy_assign --regions <whole contigs> --families run.fam.copies.tsv --copies-fa …` (`--copy-table catalog` = the legacy `gw_family_catalog` table) | `run.assign.assignments.tsv` |
| `flag` | `missing_copy_flag --scan-only` + `--from-scan` on a splice `.mmi` of `genome.ref.fa`, `--confirm truth=` the truth genome | `run.flag.missing_copy.tsv` |

Logs in `DIR/logs/` and the driver's `run.*.log`.

## `score.tsv` / `summary.md`

Long table `scenario objective metric copy value note`; `summary.md` one table per objective.

- **alignment** (reported, not judged): per copy the primary's class `own | other:<id> | absorbed_by:<id> | outside |
  unmapped`, MAPQ-0 share, AS-tied share (secondary within 0.98 of the primary AS), supplementary, soft-clip ≥ 50,
  insertion ≥ 50, junction recall of own-copy primaries, mean clipped fraction; per read in `reads.placement.tsv`.
- **o1_loci**: assembled loci ↔ reference-present copies one-to-one by shared exonic bp (scipy LSAP); per copy
  `found | split:k | merged:<ids> | missed | unexpressed`, exonic recall/precision (locus exon UNION), best transcript's
  junction P/R and `exact_chain`.
- **o1_family_denovo / o1_family_guided**: cluster per copy (members mapped to copies by span), pairwise
  sensitivity/precision/F (`lib.pairwise`), `bipartite_F` (`lib.bipartite_families`), `family_intact` (every scoreable copy of
  A in one cluster, no decoy; `n/a` with < 2 scoreable copies — de novo cannot see an unexpressed copy, neither mode sees an
  absent one), `guided_minus_denovo_F`.
- **o2**: tied reads (MAPQ 0 or AS-tied) from `assignments.tsv`: `tied_correct | tied_wrong | tied_abstain |
  tied_absent_assigned`; unique reads by the aligner's placement (the AS-tied gate skips them by design, as `sim.py tandem`
  scores them): `unique_placed_correct | _wrong | _absent`; `wrong_among_assigned`.
- **o3**: the verdict at the locus that absorbed each absent copy's reads (expected `reference_absent_candidate` there,
  `not_flagged` elsewhere); `flag_correct`; the scan's `class`, `m`, `delta`, `n_psv`, `conf_truth_identity`.

## The built-in ladder (`scenarios.py`)

22 rungs on one template, one background, 3 decoys, 50 reads per copy: `identical`, `snp_{0.001,0.005,0.01,0.02,0.05,0.1}`,
`exon_loss`, `splice_kill`, `inv_exon`, `exon_gain`, `exon_dup`, `shuffle`, `inv_intron`, `inv_whole`, `truncated`,
`conversion`, `three_copies`, `dispersed`, `unexpressed`, `absent`, `combined` (every structural op at once). Rungs that
need more exons than the template has are skipped and listed.

### 2026-10-02 run — human SNRPB (CAT chr20, 8 exons) and orangutan SPAG7 (RefSeq NC_072392.2, 7 exons)

Specs and tables: `bench/famsim/examples/ladder_{human_SNRPB,ppy_SPAG7}.json`, `.2026-10-02.tsv`. **Every rung's verify
passed in both species** (`PASS` = all claims). The pipeline columns are what the SHIPPED binaries did; no bars were
pre-registered for them (write the `PREREG_` before reading further).

| rung | human de novo / guided intact | orangutan de novo / guided intact | loci found | O2 wrong among assigned | O3 |
|---|---|---|---|---|---|
| identical | 1 / 1 | 1 / 1 | 1.0 | 0 (all A reads tied → abstain, 0.40 of reads) | — |
| snp 0.001, 0.005, 0.02, 0.05 | 1 / 1 | 1 / 1 | 1.0 | 0 | — |
| **snp 0.01, 0.1** | **0 / 0** | **0 / 0** | 1.0 | — | — |
| exon_loss, splice_kill, exon_gain, exon_dup | 1 / 1 | 1 / 1 | 1.0 | 0 | — |
| **inv_exon** | **0 / 0** | 1 / 1 | 1.0 | — | — |
| **shuffle** | **0** / 1 | 1 / 1 | 1.0 | — | — |
| inv_intron, inv_whole, truncated, conversion, three_copies, dispersed, combined | 1 / 1 | 1 / 1 | 1.0 | 0 | — |
| unexpressed | n/a / 1 | n/a / 1 | 1.0 | — | — |
| absent | n/a / n/a | n/a / n/a | 1.0 | — | `reference_absent_candidate` at A (host) in both, confirmed by the truth genome |

**What the two bold rows are (both diagnosed, neither is the simulation):**

1. **`mcl_families` on a two-member family is numerically fragile.** `snp_0.01`, `snp_0.1` (both species) and the human
   `shuffle` de novo arm all log `2 nodes, 1 edges, 0 pairs dropped` and then `0 cluster(s) >= 2 members`. The MCL
   (`annotation_families.rs::mcl`) adds self-loops of 1.0 while the edge weight is identity × coverage < 1; on two nodes the
   iteration converges toward a tie between self and mate and the convergence test (`|u − v| < 1e-7`) stops at a residual
   whose sign decides the attractor. The bit-faithful port reproduces it: `lib.mcl({("a","b"): w})` joins at w = 0.999,
   0.995, 0.98, 0.97, 0.95, 0.93, 0.85, 0.80 and SPLITS at 0.99, 0.985, 0.90, 0.75, 0.61; three nodes always join; inflation
   2.0 joins both probed weights. 57 % of genome-wide families are size 2 (`FAMILY_DEF.md`) and `--min-size 2` is a user
   decision (§6ex), so this is worth a pre-registered fix (a connected pair is a cluster; or self-loops = the column's max
   edge weight, standard MCL; or a tolerant attractor comparison). Nothing was changed here.
2. **An inverted exon inside a 2 %-divergent copy can lose the pair at `--min-shared-exon-frac 0.60`** (human SNRPB exon 5
   of 8: dropped in BOTH modes; orangutan SPAG7 exon 4 of 7: kept). The inversion splits the all-vs-all alignment into
   records and the shared-exon conjunct is evaluated on one of them.

Everything else the advisor named — fewer exons, more exons, duplicated exons, shuffled exons, inverted introns, whole-gene
inversions, partial copies, gene conversion, three copies, a copy on another contig, SNVs at 0.1–5 % — leaves the family
intact in both modes, the loci found with exact chains, zero wrong assignments, and the reference-absent copy flagged at its
host. The identical-copy rung behaves as the thesis says it must: every read of A/A′ is AS-tied and the certificate
abstains (no PSV), unique reads of the decoys are placed correctly.

## Limits and notes

- **Gorilla**: `winloci_data/GGO_genomic.gff` names RefSeq contigs (`NC_0732xx`); the only gorilla FASTAs on disk
  (`gorilla_haps/{mat,pat}.fa`) are the GenBank-named haplotype assemblies, so no gorilla ladder was run. Point the template
  at a FASTA whose contig names match the GFF and it runs unchanged. Chimp (`mPanTro3` FASTA, no GFF on disk) likewise.
- `asm20` is the pipeline's all-vs-all preset; at 10 % divergence it still aligned the copies here (identity 0.90, full
  coverage); the verifier falls back to a Python global alignment for chain identity when minimap2 finds no hit.
- Readthrough molecules, 5′-cap/polyA modelling and library-learned length distributions are out of scope (`sim.py
  chromosome rt` covers readthrough).
- Human and ape results are never pooled; run the ladder per species and report both tables.
- To add a condition: a new op in `ops.py` (+ its verify claim in `verify.py` and a test), or a new rung in
  `scenarios.rungs`.
