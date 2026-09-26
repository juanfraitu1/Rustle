# Pre-registration: the families stage writes the copy table that copy assignment consumes (parity checks)

**Written 2026-09-25 16:20, before any number below exists.** No `mcl_families --from-gtf --emit-units` output
has ever been produced (the path is untested: `src/bin/mcl_families.rs`, `--emit-units` needed `--bam`, and no
driver, bench or figure script passes `--bam` with `--from-gtf`; checked with `grep`).

## 0. What is decided, and by whom

The user decided on 2026-09-25 16:00 that the ONE default de novo family definition is the driver's `families`
stage (reads -> seeded assembly loci -> one locus representative per locus, its "positional exon sum" -> families:
`mcl_families --from-gtf`, exon-sum >= 0.60, MCL 2.8), that copy assignment consumes the SAME families, and that the
`gw_family_catalog` (`catalog` stage) becomes legacy. **That decision is not tested here and nothing below can
reverse it.** This document pre-registers only the ENGINEERING checks that the families stage now writes a copy
table in the contract `copy_assign --families/--copies-fa` and `bench/sim.py copies` read, and the descriptive
comparison with the legacy catalog. There is no accuracy claim and no decision rule on the counts.

## 1. The product

`mcl_families --from-gtf GTF --fasta G --emit-units --out P` writes, beside the unchanged `P.clusters.tsv`,
`P.loci.*` and `P.params.tsv`:

- `P.copies.tsv`: one row per member locus of every reported family (the members of `P.clusters.tsv`, in the same
  `MCL<i>` families). Columns 1-11 are the catalog's (`family_id copy_idx tid chrom start end n_exon strand n_reads
  exons max_family_identity`); further columns are appended. A copy is the locus REPRESENTATIVE (the transcript with
  the most reads, ties to the longer span: the same transcript `loci.gff3` carries): `tid` = its `transcript_id`,
  `exons` = its exons (0-based half-open), `n_reads` = its `reads`, `strand` = its strand.
  `max_family_identity` = the identity of the best families-stage edge from this locus to another member (genomic
  span alignment, `-x asm20`; NOT the catalog's exon-sum identity: a schema difference, reported).
- `P.copies.fa`: the spliced exon sum of the representative (genome bases at its exons, reverse-complemented on `-`),
  header `>family|copy|chrom:start-end|strand|nexon=N` (the catalog's).
- `P.copies.regions`: per family and contig, the copies' hull +/- 5 kb.

## 2. Development data (human A119b chr16; no held-out claim)

- de novo GTF: `/mnt/linuxdisk/tmp/rustle_figures/fig7/current/human_chr16.denovo.gtf` (the genome-wide A119b
  assembly restricted to chr16; the fig. 7 development run);
- reads: `/mnt/linuxdisk/home/juanfraitu/chr16_arm/chr16.bam`; genome: CHM13 v2.0;
- legacy catalog for comparison: `/mnt/linuxdisk/tmp/rustle_figures/fig6/chr16.cat.copies.tsv` (= `chr16_arm/on.copies.tsv`).

## 3. Checks (each PASS/FAIL; all must pass)

- **C1 byte identity.** (a) With `--from-gtf` and without `--emit-units`, every product of the new binary is
  byte-identical to the binary built from the tree before this change, same inputs. (b) With `--emit-units`,
  `clusters.tsv`, `loci.tsv`, `loci.gff3`, `loci.fa`, `loci.paf` are byte-identical to (a); `params.tsv` differs only
  in the `emit_units` row and in rows appended after the last existing row. (c) The annotation mode (`--paf --gff`,
  with and without `--bam` read-chain units) is byte-identical to the old binary on the same inputs.
- **C2 schema.** Header columns 1-11 equal the legacy catalog's header, in order; every row passes
  `catalog_input::parse_copies_tsv` rules (exons ascending, disjoint, reconstruct start/end and `n_exon`; strand
  `+`/`-`; locus extent contains the span); every FASTA record's header matches its row and its sequence equals the
  genome at the row's exons (reverse-complemented on `-`), checked for EVERY copy with `samtools faidx`.
- **C3 copy_assign loads and completes.** `copy_assign --families P.copies.tsv --copies-fa P.copies.fa` on the real
  chr16 BAM, restricted to a small subset of families (at most 10, including the largest), exits 0.
- **C4 simulation.** `bench/sim.py copies` runs on the families table (same subset) against a chr16-only splice
  index, and `copy_assign --families` on the simulated BAM exits 0.

## 4. Descriptive only (reported, no bar, no decision)

Copies and multi-copy families: families table vs legacy catalog; copies of one table whose exons overlap a copy of
the other (both directions); strand and exon-count distributions.

## Amendments

(none)
