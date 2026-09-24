# Pre-registration — is our per-transcript quantification comparable to StringTie / FLAIR / isoseq collapse?

**Written 2026-09-23 (§6zc), before any correlation is computed.** User: *"ensure the quantification is
comparable to StringTie's or the other tools for the matching transcripts that we also find … in TPMs and
RPKMs if needed."* Item 1 of `docs/PENDING_2026-09-23.md`.

## Objects

- **Ours:** the new-default genome-wide `--assemble-only` output (human A119b, 25 contigs; gorilla, 26
  contigs), `reads "N"` per transcript = distinct-coordinate primary molecules carrying the exact chain
  (de-duplicated, r1058). Also the same run with `--keep-coordinate-duplicates` (every primary molecule),
  because the tools do not de-duplicate.
- **Tools (lab runs, same BAMs):** StringTie 3.0.1 `cov` / `FPKM` / `TPM` attributes; FLAIR 3.0.1
  `*.isoform.counts.txt` (reads assigned by its quantify step; key = `transcript_id + "_" + gene_id`, r1060
  note); isoseq collapse `*.collapsed.abundance.txt` `count_fl`.
- **Join:** strand-aware multi-exon intron chain, exact (the r1065 key). Single-exon transcripts are excluded
  (no chain key; their ends are not comparable across tools).

## Quantities

For each shared chain: our `reads`; our count-based TPM = reads / Σ reads × 10⁶ over our whole output; our
RPKM = reads / (spliced kb × Σ reads / 10⁶) — reported because it was asked for, with the standing caveat
that a long-read count is already a molecule count and length-normalising it is wrong for this data (§6r5:
count-based ρ 0.879 vs StringTie, length-normalised 0.714 on chr20). Tool side: isoseq `count_fl`, FLAIR
count, StringTie `cov` (reads per bp, its molecule-count proxy) and StringTie `TPM`.

## Metrics — committed now

Per species, per tool, on the shared chains: **Spearman ρ** (primary) and Pearson r of log1p values, the
number of shared chains, the median ratio ours / tool, and the same restricted to the multi-copy families the
thesis is about (transcripts whose chain is annotation-exact at a gene named `NPIP*`, `TBC1D3*`, `NBPF*`,
`GOLGA*`, `LOC*` gorilla equivalents are reported as a class).

| outcome | verdict |
|---|---|
| ρ ≥ 0.85 vs StringTie `cov` AND vs isoseq `count_fl`, on ≥ 1,000 shared chains, both species | ⭐ **COMPARABLE** |
| ρ ≥ 0.75 on both, or ≥ 0.85 on one species only | ⚠ **CLOSE, with a stated systematic** |
| below | ⛔ **NOT COMPARABLE** — find the cause before showing the advisor |

**Predicted, before looking:** ⭐ against isoseq (both are molecule counts of the same reads; the de-duplicated
count sits below `count_fl` at deep loci by the 25-39% duplicate rate, a monotone effect that Spearman
ignores), ρ ≈ 0.85-0.90 against StringTie `cov` (its coverage model spreads reads across isoforms), and lower
against FLAIR (its counts include partial-alignment assignments, r1057). RPKM will correlate worse with
everything, as measured before.

I will not change the join, the quantities, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-23) — ⚠ **CLOSE, with one stated systematic that a flag removes**

Shared strand-aware multi-exon chains, genome-wide, tool counts SUMMED over a tool's isoforms that share the
chain (isoseq keeps several PB ids per chain), Spearman ρ / Pearson r on log1p (`/mnt/linuxdisk/tmp/gw22/quant/parity.py`):

| species | our counts | vs isoseq `count_fl` (n) | vs StringTie `cov` (n) | our TPM vs StringTie TPM | vs FLAIR count (n) |
|---|---|---|---|---|---|
| gorilla | de-duplicated (default) | 0.786 / 0.78 (57,544) | **0.898** / 0.94 (43,079) | 0.898 / 0.94 | 0.786 / 0.87 (49,853) |
| gorilla | all molecules (`--keep-coordinate-duplicates`) | 0.786 / 0.77 (50,757) | **0.906** / 0.95 (41,337) | 0.906 / 0.96 | 0.824 / 0.91 (43,066) |
| human | de-duplicated (default) | 0.803 / 0.85 (163,193) | 0.832 / 0.89 (86,058) | 0.832 / 0.90 | 0.717 / 0.81 (119,389) |
| human | all molecules | **0.833** / 0.88 (175,191) | **0.836** / 0.91 (97,762) | 0.836 / 0.93 | 0.699 / 0.83 (121,577) |

**The systematic is the de-duplication (r1058).** With the default counts the median ratio ours / isoseq
falls with depth — 1.00 at ≤ 4 FL reads, 0.75 at 16, 0.33 at ≥ 256 (human) — because coordinate-identical
molecules are counted once. With `--keep-coordinate-duplicates` the ratio is **1.00 at every depth bin on
both species** (0.98 at ≥ 256 on human): our molecule count IS isoseq's FL count. StringTie's `cov` sits
below ours by a constant (median ours / cov 0.87 gorilla, 0.80 human) with no depth trend; our count-based
TPM tracks its TPM at 0.91-0.96 Pearson.

By the bar: StringTie clears 0.85 on gorilla (0.906) and misses it on human by 0.014 (0.836); isoseq sits at
0.79 / 0.83 — **⚠ CLOSE**, ρ ≥ 0.75 everywhere. The residual is at the low-count end (most shared chains carry
2-4 reads, where a one-read difference in how a tool assigns reads among near-identical isoforms moves rank
a lot) — but the ≥ 10-read restriction below shows that is only part of it. RPKM, as predicted, correlates worse with everything and is
the wrong normalisation for molecule counts. Multi-copy family chains (NPIP/TBC1D3/NBPF/GOLGA, n 52-71 on
human): ρ 0.82-0.93 vs every tool.

**For the advisor:** quantify with `--keep-coordinate-duplicates` (the tools' convention) and compare
count-based TPM; the de-duplicated default is the better support criterion for assembly (r1066) but not
the comparable quantity.

**Restricted to chains the tool counts at ≥ 10 reads (all-molecule counts):** gorilla isoseq 0.793 (n 14,560),
FLAIR 0.804, StringTie `cov` 0.880, TPM 0.890; human isoseq **0.852** (n 41,659), FLAIR 0.656, StringTie
`cov` 0.799, TPM 0.802. Depth lifts only the isoseq/human cell over the bar; StringTie's cells FALL — the
residual against StringTie is its coverage model dividing reads among overlapping isoforms differently from
an exact-chain molecule count, not count noise. Both are "comparable"; they are not the same quantity, and
the write-up says so rather than tuning toward either.
