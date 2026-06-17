# Two synthetic benchmark datasets (canonical "ideal cases" for advisor meetings)

Both use a deterministic full-length HiFi transcript simulator (`bench/sim_reads.py`): each read = one
transcript + HiFi errors (substitutions + rare short indels), a KNOWN error rate (what the
identifiability theorem needs). Large artifacts live in `/home/juanfra/winloci_scratch/{topup,sim5x}/`
(not in the repo); the repo holds the builders, figures, and summary tables.

---

## Dataset 1 — "ideal coverage" GGO (top up only what the real IsoSeq lacks)
`bench/build_topup.py` → `/home/juanfra/winloci_scratch/topup/GGO_ideal.bam`

The real GGO IsoSeq is coverage-uneven (median 43 reads/gene; **2,223 genes with ZERO reads**), and
coverage limited every prior analysis (headroom = 0 was partly coverage-bound; both-copies-expressed
was only ~3%). This dataset removes coverage as a confound by **simulating only what's lacking**: every
gene under 40× is topped up to 40× from its representative transcript, mapped to the genome
(`minimap2 -ax splice:hq`), and merged with the real BAM.

- **11,225** under-covered genes topped up; **344,868** synthetic reads; **100%** map back.
- **`GGO_ideal.bam` = 1.8 GB** (real + synthetic). 96% of topped-up genes now reach ≥40×
  (residual = long genes where truncated reads don't span the full window).
- Ground truth: synthetic read names `SIMTOPUP|<gene>|<i>` + `topup_truth.tsv`.
- Disk-light: simulate → map → delete the 1.7 GB FASTQ/SAM; keep only the BAMs.

**Use:** re-run the recall / allele-specific / copy analyses on `GGO_ideal.bam` to see what's findable
when coverage is *not* the limiting factor — the "ideal" upper bound to compare the real data against.

![coverage](topup_coverage.png)

---

## Dataset 2 — "5 equally good places" (the copy-assignment identifiability benchmark)
`bench/build_sim5x.py` → `/home/juanfra/winloci_scratch/sim5x/` (ref.fa + bam per K + summaries)

A real gene (**AASDHPPT**, 6 exons, 2,947 bp mRNA) replicated into **5 near-identical copies in tandem
in the reference**, so a read has *five equally good places to put it* — the regime the advisor asks
about. A **PSV divergence ladder** (K = 0 identical … 8 private exonic PSVs per copy) maps the
identifiability boundary; spliced HiFi reads (40/copy, ground-truth copy in the name) are mapped back
with `minimap2 -ax splice:hq`.

**The result — coordinates cannot assign; PSVs can, iff ≥K columns clear the error floor:**

| K (PSVs/copy) | % reads MAPQ-0 (coords can't assign) | copies identifiable /5 | PSV-assignment accuracy |
|---|---|---|---|
| 0 (identical) | **100%** | 0 | 0.20 (random) |
| 1 | 40% | 3 | 0.80 |
| 2 | 0% | **5** | **1.00** |
| 4 | 0% | 5 | 1.00 |
| 8 | 0% | 5 | 0.99 |

- **K=0 (identical copies)** is the hard extreme: minimap2 leaves 100% of reads at MAPQ 0, and there is
  *no PSV information* to recover from — assignment is **information-theoretically impossible** (accuracy
  = random 1/5). This is the canonical case to point at.
- **K=1**: 4 bases can't separate 5 copies (copies collide) → only 3/5 identifiable.
- **K≥2**: 5 copies become fully separable (4² = 16 ≥ 5) → both minimap2 and PSV-assignment resolve them.
  The boundary is at **⌈log₄ N⌉ = 2** PSV columns for N=5 copies.

**The error floor (2nd theorem axis, at K=4):** PSV-assignment accuracy 1.00 at HiFi error (e=0.003) →
0.985 (1%) → 0.905 (10%) → 0.805 (15%). PSVs assign copies only while they clear the per-base error.

![identifiability](sim5x_benchmark.png)

**Use:** the controlled testbed for copy-assignment / `copy_split` / the identifiability theorem — vary
the number of copies, PSV count, coverage, and error to show exactly when "which copy did this read
come from?" is answerable and when it is not.

## Reproduce
- `python3 bench/gene_coverage.py` (per-gene real coverage → defines the top-up set)
- `MINIFORGE python bench/build_topup.py` (Dataset 1) ; `python3 bench/topup_fig.py`
- `MINIFORGE python bench/build_sim5x.py` (Dataset 2) ; `python3 bench/sim5x_fig.py`
