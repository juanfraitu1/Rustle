# Re-running the RefSeq benchmarks on CAT/Liftoff v2.0: protocol (2026-10-01)

Written before any benchmark was re-run. Follows `docs/ANNOTATION_CAT_DEFAULT_2026-10-01.md` (user decision: CAT/Liftoff v2.0 is
the default human annotation; registered RefSeq results stay as registered; each benchmark is re-run on CAT and both are reported).
Annotation: `winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz` (gene `Name=` = CAT gene id; `gene_name=` = symbol).
Gorilla is untouched. Nothing in the canonical repo is edited; changed copies of scripts live on branch `machine2/soto-evidence`
(repo scripts) or next to the original scratch instruments with a `_cat` suffix (scratch scripts), originals untouched.

## Rulings (fixed now)

- **R1 Re-keying a RefSeq-defined truth copy.** A copy maps to the CAT gene sharing the most exonic bases with it on the same strand
  (`chm13v2.0_CAT_Liftoff.refseq_map.tsv`; ties: larger Jaccard). Never by name. A copy with no CAT gene ("none") is dropped from the
  CAT truth and listed. Other CAT records overlapping the copy (nested lncRNAs, read-through models, Liftoff duplicates) are not
  copies; they are listed. A RefSeq record without exon lines (NPIPB14P, TBC1D3P1/P3/P4/P7) maps by its gene span to the CAT gene with
  the largest span overlap on the same strand, and becomes scorable if that gene has exons.
- **R2 Members defined by RefSeq text.** Where a benchmark selected members by the RefSeq `description` (e.g. "TBC1 domain family member
  3", "nuclear pore complex interacting protein"), the CAT member set is the R1 image of the RefSeq member set; the CAT genes whose
  `gene_name` starts with the family symbol and contains no `-` are reported beside it, never silently substituted.
- **R3 Read-throughs.** CAT has no `description` field. A CAT gene is a read-through when its `gene_name` is `A-B` with both `A` and `B`
  CAT gene names, or when it is the R1 image of a RefSeq read-through.
- **R4 HGNC truth.** Via CAT `source_gene` (Ensembl gene id) -> HGNC `ensembl_gene_id` -> `gene_group_id`; genes without a match are
  listed.
- **R5 Soto joins** use S1C `Gene ID` (CHM13_G / LOFF_G, = CAT gene id), never `Gene Name`.
- **R6 Reporting.** Each benchmark reports RefSeq (registered) and CAT side by side with the same scorer; denominators that change
  (25 vs 26 NPIP copies; TBC1D3 pseudogenes that gain exons) are stated. No threshold, rule or parameter is re-tuned on CAT.

## Order (dependencies)

1. Re-key the NPIP / TBC1D3 truth tables (Dishuck 27-row table; TBC1D3 members).
2. Copy recovery (scoring only; tools ran annotation-free).
3. Held-out families (PREREG 2026-09-20): CAT regions, all-vs-all, `mcl_families`, `score.py heldout` joined by Gene ID.
4. SQANTI3: CAT reference GTF, same six chromosomes and arms.
5. Layer-order / nested-lattice `light/` tables rebuilt from CAT (on a copy of ROOT).
6. O1 L-level scoring (lattice engine, family certificates) on the rebuilt tables.

## Known before starting (from the mapping, 2026-10-01)

RefSeq NPIPB3 (chr16:21,337,400-21,360,419) has no CAT gene: the CAT NPIP truth has 25 of 26 chr16 copies (+ NPIPB1P chr18). CAT names
are permuted relative to RefSeq (CAT NPIPB3 = RefSeq NPIPB5; TBC1D3 D/E/K rotated). RefSeq NPIP / TBC1D3 copies vs CAT: 18 strong, 15
partial, 1 weak, 7 none (before R1's span rule for exon-less records).
