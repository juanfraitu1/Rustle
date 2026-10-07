# Ideal few-copy families for O1: input-side criteria (human, CAT/Liftoff v2.0, A119b) — FROZEN before any count

Written 2026-10-06 before any criterion value, filter count or shortlist was computed. Human only. Nothing in this selection reads a pipeline
product: no `*.fam.clusters.tsv`, no `family_score` output, no `outcomes.tsv`, no assembled GTF, no AS table. Inputs: the Compara family truth,
`refseq_map.tsv`, the CAT/Liftoff `genes.tsv` and slim GFF3, the CHM13 v2.0 FASTA, and region queries of `A119b.t2t.bam`.

What was looked at before freezing (no per-family value): the family-size histogram of the truth file (377 of 426 families have 2-4 listed genes: 292 / 60 / 25),
the vocabulary of `refseq_map.tsv` quality (strong / partial / weak / none) and of the CAT transcript biotypes, the BAM header (minimap2 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`,
tags NM ms AS nn ts tp cm s1 de rl), `docs/PREREG_ideal_expression_2026-10-06.md` and `bench/ideal_expression/strata.py` (the E0 rule), and two format probes of one NPIP region (chr16:14938000-14943000, record tags only).

## Universe

Compara Primates families of `/mnt/linuxdisk/tmp/rustle_figures/families_gw/species/human/compara.Primates.families.tsv` (columns Gene Name, Family ID, Contig) with exactly 2, 3 or 4 listed genes (N0 = 377).
The file lists only Compara genes that were matched to a CHM13 RefSeq protein-coding gene, so a listed size of 2-4 can hide unlisted members (a limitation, not a criterion).

## Member mapping (feeds C1; an unmapped member is recorded with its reason)

For member (symbol S, contig K): RefSeq candidates = rows of `chm13v2.0_CAT_Liftoff.refseq_map.tsv` with `refseq_name == S`, `chrom == K`, `refseq_biotype == protein_coding`.
0 rows -> UNMAPPED:symbol_not_found. More than 1 row -> UNMAPPED:ambiguous_symbol (never pick one). 1 row but empty `cat_gene` (quality none) -> UNMAPPED:no_cat_gene.
Otherwise the member's CAT/Liftoff gene g = `cat_gene` of that row (the gene sharing the most exonic bp on the same strand, `bench/annotation/cat_setup.py`; names are never used to match).

## Reference transcript (RT)

Among the transcripts of g in the slim GFF3 with `transcript_biotype=protein_coding`, the one with the largest spliced length (sum of merged exon lengths); ties: more exons, then the smaller transcript ID.
Exons of a transcript are sorted and overlapping or adjacent ones merged; an intron is (end of upstream exon, start of downstream exon), 0-based half-open; its motif is read from `chm13v2.0.fa` on the gene strand
(reverse-complement for `-`, upper-cased): CANONICAL = GT-AG, GC-AG or AT-AC (donor-acceptor of the first two and last two intron bases).

## Criteria (a family passes a criterion only if EVERY member / pair passes; NA counts as fail)

- **C1 annotation:** every member is mapped; every member's `refseq_map` quality is `strong` (shared exonic bp >= half of BOTH exon unions); every member's CAT/Liftoff gene has `gene_biotype == protein_coding`;
  the members' CAT/Liftoff gene IDs are pairwise distinct; every member has an RT; the RT has >= 3 exons; every intron of the RT is >= 50 bp and CANONICAL.
- **C2 not entangled:** for every member's gene g (exon union over all its transcripts) and every OTHER CAT/Liftoff gene h of `genes.tsv` on either strand (any biotype, any source, other family members included),
  the shared exonic bp |union(g) ∩ union(h)| is < 100.
- **C3 separate loci:** for every pair of members on the same chromosome, the gap between their CAT/Liftoff gene spans (`start0`, `end` of `genes.tsv`) is >= 1000 bp (gap = max(start) - min(end); an overlap is negative). Members on different chromosomes pass.
- **C4 homologous, not identical:** RT spliced sequences (transcript orientation, from the FASTA) of every pair aligned with `minimap2 -x asm20 --cs -t 1` (target = the longer sequence, query = the shorter; ties by CAT gene ID), PAF records with `tp:A:P` or `tp:A:I` only.
  cov = length of the union of query intervals / length of the shorter sequence; identity = sum(PAF col 10 matches) / sum(PAF col 11 block length). PASS iff, for EVERY pair, cov >= 0.50 AND 0.90 <= identity < 0.999
  (0.90 = the usual segmental-duplication floor and inside asm20's design range; 0.999 = about one difference per kb, the point where HiFi error (about 1e-3) can no longer tell the copies apart).
  A pair with no alignment has cov 0 and fails. Descriptive tier of the family's minimum identity (never a gate): close >= 0.98, mid 0.95-0.98, far 0.90-0.95.
- **C5 expressed in A119b (E0 rule):** region queries only, `pysam fetch(contig, start0, end)` over each member's CAT/Liftoff gene span, every alignment (primary, secondary, supplementary; any MAPQ; unmapped skipped).
  Valid chains of a member = the intron tuples of those transcripts of g (any biotype) that have >= 2 introns after merging and whose every intron is >= 50 bp and CANONICAL.
  A read alignment carries chain c iff the list of its CIGAR `N` junctions of length >= 50 bp lying wholly inside the member's span equals c exactly (the E0 convention of `strata.py`: `M = X D` advance the reference, `N` is an intron, `I S H` do not).
  Pools, counted as DISTINCT READ NAMES per member and chain: **P1** = the alignment is primary (flag & 2308 == 0); **P2** = P1 plus secondary alignments (flag 256) whose read's primary alignment was found among the alignments fetched
  for the family's member spans and whose AS >= 0.98 x the best AS of that read (max AS over its non-supplementary alignments found there) - the pipeline's seeding pool; **P2u** = secondary alignments carrying the chain whose read's primary alignment was
  NOT found in the family's spans (primary elsewhere, status unknown: reported, never counted in P2, never gating); **P3** = any alignment, supplementary included.
  A member's support n = max over its valid chains of |P2 reads of that chain| (the support chain: ties by more introns, then lexicographic). **PASS iff n >= 3 for EVERY member.** A member without a valid chain has n = 0.
  Paralog-tie fraction of a member = (P2 reads of its support chain whose alignment at the member is NOT a primary with MAPQ > 0, i.e. secondary or MAPQ 0) / n; also printed: the share of those whose primary alignment lies on a different member of the family.
  Family tie fraction = pooled over members (sum tied / sum n).
  Descriptive, never gating: n_P2 with 5 bp junction tolerance (same number of inside junctions, every junction end within 5 bp of the chain) to expose annotation-boundary mismatches.

## Overall pass, order of the counted steps, ranking

PASS-ALL = C1 and C2 and C3 and C4 and C5. Filter counts are reported cumulatively in the order C1 -> C2 -> C3 -> C4 -> C5 from N0, and as independent per-criterion pass counts among the families on which the criterion is evaluable (C1: all families; C2, C3, C5: every member has a CAT/Liftoff gene; C4: every member has an RT).
Shortlist = PASS-ALL families ranked by (1) minimum over members of the P2 support n, descending; (2) minimum over members of the P1 count of the support chain, descending; (3) fewer members; (4) Family ID ascending.
Near-miss list (descriptive): families that fail exactly one criterion, with which one. A rule above is never edited after the first count; a code bug fix that changes a value is logged in `amendments.md` with the reason and the before/after.

## Outputs (all under the scratch directory)

`criteria.md` (this file) + `criteria.sha256` (its hash, taken before the first count), `inputs.tsv` (one row per family: every value and pass/fail per criterion, rank), `members.tsv` (one row per member), `pairs.tsv` (one row per member pair),
`unmapped.tsv`, `shortlist.tsv`, `near_miss.tsv`, `counts.json`, `provenance.txt` (input sizes, mtimes, sha1 of small inputs and of the scripts), the scripts and the minimap2 FASTA/PAF files.

## Declared limitations of the design (not gates)

(1) the family size is the number of listed genes, not the number of Compara members; (2) there is no genome-wide homology search (unannotated or non-coding paralogs elsewhere are invisible); (3) P2 resolves a secondary alignment's primary only inside the family's own spans,
so P2 can be a lower bound (P2u shows how much); (4) exact-chain matching ignores reads with 5' truncation or boundary jitter (the tolerance column shows how much); (5) one testis library (A119b), human only; (6) annotation-driven chains: CAT/Liftoff projections can be wrong;
(7) selecting on input-side ideality builds a positive-heavy set, so an O1 score on it is a mechanism check, not a recall estimate, and the shortlist should be split into dev and held-out before any O1 output is read.
