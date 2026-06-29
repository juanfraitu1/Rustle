# Facultative long-read phasing (`copy_assign --phase`) — dependency-free

The PSV harness already phases reads into copy-haplotypes internally (`copy_split::split_readchain_by_psv`
+ `psv_linkage::assign_read_to_copy`). `--phase` surfaces that as a first-class, **opt-in,
dependency-free** phasing output — no external phaser (WhatsHap/HapCUT2), no neural variant caller
(DeepVariant/Clair3), no CNN. It is the N-copy generalisation of read-backed phasing: a long molecule
spans multiple PSVs, the linked pattern assigns it to a haplotype, and the formal object is the
minimum path-cover of the PSV graph (MCC = χ(H), the identifiability theorem).

## Outputs (written only under `--phase`; all other outputs unchanged → additive/default-off)

- **`<out>.phase.gfa`** — the phasing AS A VARIATION GRAPH (self-contained GFA1, loadable in
  Bandage / `vg`): each PSV column is a **bubble** (one `S` segment per allele, carrying the base +
  `PO:i:<genomic_pos>`); each copy is a **path** (`P` line) threading the bubbles. Copies that share an
  allele share the node (a bubble anchor) and fork where they differ — so the graph topology *is* the
  copy structure. PSV-identical copies collapse to **one path** = the K-frontier rendered as graph.
- **`<out>.phase_blocks.tsv`** — one **phase set (PS)** per family:
  `block_id  chrom  n_haplotypes  n_psv_sites  n_reads_phased  n_unphased`
- **`<out>.phased_haplotypes.tsv`** — each haplotype's phased alleles:
  `block_id  haplotype  copy_tid  n_support_reads  variants` where `variants` = `pos:allele;pos:allele;…`
- **`<out>.phased_reads.tsv`** — the **haplotag (HP)** per read:
  `read_name  block_id  haplotype  n_psv_spanned  margin  status`
  (`haplotype = -1` when the read is Ambiguous/Tied = unphaseable = the K-frontier).

`block ≙ PS`, `haplotype ≙ HP`, `bubble ≙ PSV`, `path ≙ copy` — standard phasing/VG terms. The `.gfa` +
`.phased_reads.tsv` are jointly self-contained: the graph holds copies-as-paths-through-PSV-bubbles, the
haplotag maps each read to its copy/path.

### GFA demo (NC_073242.2 family: 3 copies, 43 PSV columns)

```
phase.gfa: 86 bubble-nodes, 84 links, 3 copy-paths
S  CAFAM0_c0_C  C  PO:i:3790785     S  CAFAM0_c0_G  G  PO:i:3790785   <- the bubble (2 alleles)
P  CAFAM0_copy0  c0_C+,c1_G+,...  ┐ identical path (PSV-identical copies collapse = K-frontier)
P  CAFAM0_copy1  c0_C+,c1_G+,...  ┘
P  CAFAM0_copy2  c0_G+,c1_C+,...    the divergent copy forks at every bubble
```

## Demo (NC_073242.2:3771193-3799186, a clean tandem family)

```
phase_blocks:  CAFAM0  NC_073242.2  3 haplotypes  43 PSV sites  18 phased  0 unphased
haplotypes:    hap2  2 reads  3790785:G;3790787:C;…   (PSV-distinguishable; margin 330)
               hap0  8 reads  3790785:C;3790787:G;…   ┐ identical PSV alleles, resolved by
               hap1  8 reads  3790785:C;3790787:G;…   ┘ LOCUS/read-coherence axis (margin 20)
phased_reads:  m64076…/161284984/ccs  CAFAM0  hap0  45 spanned  margin 20.0  assigned
```

Two axes at once: PSV linkage (hap2) **and** read-coherence/locus (hap0 vs hap1, PSV-identical copies).
On a pathological size-heterogeneous family (PCDHB container-locus, 116 kb backbone + ~4 kb paralogs) it
phases **0/15,386** reads — correctly **abstaining** rather than guessing (DAZ3/Canzar discipline).

## Properties

- **Dependency-free**: pure Rust over data already computed (`copy_psv_alleles`, `psv_col_pos`,
  `assignments`); no new crates, no external tools, no ML (grep-confirmed: zero torch/onnx/sklearn/CNN in
  the harness).
- **Opt-in / additive**: gated behind `--phase`; without it, the families/assignments/quant outputs are
  byte-identical.
- **Confidence-gated**: a read is phased iff it clears the decisive-margin gate τ = ln((1−p)/p)
  (default 6.9 ≈ p 1e-3, the Eichler AS≥10 analog); unphaseable reads are emitted as `-1`, never guessed.
- **Phases copies, not alleles**: it phases into paralog copy-haplotypes (the thesis problem), not the
  two parental allelic haplotypes; allele-level het is the separate ASJ machinery.

Implementation: `src/bin/copy_assign.rs` (`--phase` flag + emission). Build: `cargo build --release --bin copy_assign`.
