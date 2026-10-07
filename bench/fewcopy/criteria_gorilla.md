# Ideal few-copy families for O1: input-side criteria, GORILLA (RefSeq GCF_029281585.2-RS_2024_02 / OR6737 testis IsoSeq) — FROZEN before any count

Written 2026-10-07 (local 00:0x PDT) before any criterion value, filter count or shortlist was computed on gorilla. Adapted from the frozen human criteria `bench/fewcopy/criteria.md`
(sha256 `89363117a885718695d70c6d71a78b75f88c31a62421904c505aeba34c7c1571`, re-checked identical by `sha256sum` today). Gorilla only: human numbers are never pooled with these.
Nothing in this selection reads a pipeline product: no `*.fam.clusters.tsv`, no `family_score` output, no assembled GTF, no AS table, no StringTie / FLAIR / isoseq output, no register text on gorilla outcomes.
Inputs: the projected truth `gorilla_truth.tsv` (sha256 `9e6cb80ad518f5282da3b23bd4a9f4cd89ff1bfbd39f9160b6c670d8309eb911`) and `gorilla_truth_families.tsv`, the RefSeq GFF `GGO_genomic.gff`, the FASTA `GGO.fasta` (+ `.fai`),
and region queries of `GGO_mm.bam` (+ `.bai`). Paths: truth = SCRATCH/gorilla_truth*.tsv; GFF = /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff; FASTA = /mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta;
BAM = /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam (SCRATCH = /tmp/claude-1000/-mnt-c-Users-jfris-Desktop/cde4361b-b077-41a7-8dcf-9487f51a91c7/scratchpad/fewcopy_gorilla/).

What was looked at before freezing (no per-family value): the human criteria, the three required docs, `feedback_metric_traps.md`, the human instrument scripts (`fewcopy_inputs.py`, `verify_independent.py`, `cross_links.py`, `amendments.md`, `README_inputs_side.txt`, `provenance.txt`, `checks.txt`, `counts.json`
of /mnt/linuxdisk/tmp/fewcopy_2026-10-06/; no outcomes file); the header and first rows (CF1, CF2) of the truth tables and the aggregate truth summary given in the task; the truth agent's `notes_progress.md` (exposure-ledger notes, no family values);
the GFF header, its genome-wide feature-type counts (gene 34,114; pseudogene 7,079; mRNA 80,997; lnc_RNA 10,297; transcript 3,725; exon 1,206,330; CDS 1,022,353; ...), gene_biotype vocabulary (gene: protein_coding 22,650, lncRNA 7,330, snRNA, snoRNA, tRNA, miRNA, rRNA, V_segment, ncRNA, misc_RNA, C_region;
pseudogene: pseudogene 7,027, transcribed_pseudogene 52), the exon-parent vocabulary (every exon has exactly one Parent; Parent is an `rna-*` or `id-*` transcript ID, or a `gene-*` ID for 14,597 direct exons such as pseudogenes; every transcript-level feature has Parent `gene-*`: a strict two-level hierarchy);
and the record tags of 4 alignments of one arbitrary region of the BAM (NC_073224.2:50,000,000-50,300,000; tags NM ms AS nn ts tp cm s1 s2 de rl; `--eqx` CIGARs with `=` `X` `N`); `.fai` and BAM @SQ identical (truth agent, `truth_work/verify/`). Whether that BAM region overlaps a family member was not checked.

## Data-forced changes relative to the human criteria (everything else is identical; EVERY threshold is kept: 50 bp, 100 bp, 1000 bp, 0.50, [0.90, 0.999), 0.98, 3 reads, 5 bp)

1. Annotation = the single RefSeq gorilla annotation. There is no CAT/Liftoff mapping step and no `refseq_map` quality tier. A member is the RefSeq gene-level record `gene` with `gene_biotype=protein_coding` whose `Name` equals the human symbol (exact symbol only, from `gorilla_truth.tsv`, `member_status == 1to1`).
   The human clauses "refseq_map quality is strong" and "CAT gene has gene_biotype protein_coding" are replaced by "the member is 1to1 projected" (the biotype clause then holds by construction).
2. RT = the longest transcript among the member gene's child features of type `mRNA` (instead of `transcript_biotype=protein_coding`); the tie-break is unchanged.
3. Gene IDs = the GFF gene ID (`ID=gene-...`) and the NCBI GeneID (Dbxref); there are no CAT gene IDs.
4. Spans and exon unions come from the GFF exon features through the Parent chain, never from gene-record coordinates and never joined on gene NAME (the human Amendment 1 lesson: span = extent of the exon union).
5. Valid chains use every transcript of the gene whatever its feature type (the human "any biotype").
6. The GFF is read by a single pass (no tabix); contigs are the NC_ names, identical in the GFF, the FASTA and the BAM.
7. Family IDs (CFn) are the human Compara family IDs carried by the truth table; they are labels only, never a join key to human results.

## Universe

N0 = the 377 families of `gorilla_truth_families.tsv` (Compara Primates families with 2, 3 or 4 listed human genes; 864 members; sizes 292 / 60 / 25). Family status in the truth (symbol projection): whole = every member 1to1, partial = some but not every member 1to1, none = no member 1to1.
Inputs are computed for the whole families only (member values for C2-C5 exist only where every member has a gene). The partial and none families fail C1 (a member is absent) and are listed separately in `partial_families.tsv`, with each member's state (P projected, L absent by symbol but LOC-named candidates by name evidence, N absent without evidence). LOC evidence never decides membership.

## Gene model objects (computed from the GFF)

- Gene-level records: features of type `gene` or `pseudogene` (all have IDs `gene-*`). A transcript of gene G = a feature whose Parent is G and that owns at least one exon (exon.Parent == its ID), whatever its type (mRNA, lnc_RNA, transcript, V_gene_segment, ...). Exons whose Parent is the gene itself form one direct-exon transcript of that gene (never an RT).
- Exons of a transcript are sorted and overlapping or adjacent ones merged; an intron is (end of upstream exon, start of downstream exon), 0-based half-open. Its motif is read from `GGO.fasta` on the gene strand (reverse-complement for `-`, upper-cased): CANONICAL = GT-AG, GC-AG or AT-AC (donor-acceptor of the first two and last two intron bases). An intron shorter than 4 bp has no motif and is not canonical.
- Exon union of a gene = the merge of all exons of all its transcripts. **Span of a gene = (first start, last end) of its exon union** (0-based half-open). Used for C2 candidates, C3 and C5.

## Member mapping (feeds C1)

member (family, human symbol) -> `member_status` of `gorilla_truth.tsv`. 1to1 -> gene record = `gorilla_gff_id`, which must exist in the GFF as a `gene` with `gene_biotype=protein_coding` (else UNMAPPED:record_missing). absent (L or N) -> UNPROJECTABLE (reason = the truth flags). multi / nonprotein_only -> UNPROJECTABLE (none exist in the truth, 0 and 0).

## Reference transcript (RT)

Among the `mRNA` transcripts of the member gene, the one with the largest spliced length (sum of merged exon lengths); ties: more exons, then the smaller transcript ID (string order of the GFF ID). A member with no `mRNA` has no RT.

## Criteria (a family passes a criterion only if EVERY member / pair passes; NA counts as fail)

- **C1 annotation:** every member is 1to1 projected and its record exists; the members' gene IDs (GFF gene ID and GeneID) are pairwise distinct; every member has an RT; the RT has >= 3 exons (after merging); every intron of the RT is >= 50 bp and CANONICAL.
  Sub-clauses are also counted separately (C1a whole projection, C1b distinct IDs, C1c has RT, C1d RT >= 3 exons, C1e RT introns).
- **C2 not entangled:** for every member's gene g (exon union over all its transcripts) and every OTHER gene-level record h of the GFF (type gene or pseudogene, any biotype, either strand, other family members included, on the same contig, candidates chosen by exon-union extent),
  the shared exonic bp |union(g) ∩ union(h)| is < 100.
- **C3 separate loci:** for every pair of members on the same contig, the gap between their spans is >= 1000 bp (gap = max(start0) - min(end); an overlap is negative). Members on different contigs pass.
- **C4 homologous, not identical:** RT spliced sequences (transcript orientation, from `GGO.fasta`) of every pair aligned with `minimap2 -x asm20 --cs -t 1` on two kilobase FASTA files (target = the longer sequence, query = the shorter; for equal lengths the member with the smaller GFF gene ID is the target), PAF records with `tp:A:P` or `tp:A:I` only.
  cov = length of the union of query intervals / length of the shorter sequence; identity = sum(PAF col 10 matches) / sum(PAF col 11 block length). PASS iff, for EVERY pair, cov >= 0.50 AND 0.90 <= identity < 0.999. A pair with no alignment has cov 0 and fails.
  Descriptive tier of the family's minimum identity (never a gate): close >= 0.98, mid 0.95-0.98, far 0.90-0.95. The genome index is never used.
- **C5 expressed in OR6737 (E0 rule):** region queries only, `pysam fetch(contig, span0, span1)` over each member's span, every alignment (primary, secondary, supplementary; any MAPQ; unmapped skipped).
  Valid chains of a member = the intron tuples of those transcripts of g (any feature type) that have >= 2 introns after merging and whose every intron is >= 50 bp and CANONICAL.
  A read alignment carries chain c iff the list of its CIGAR `N` junctions of length >= 50 bp lying wholly inside the member's span equals c exactly (`M = X D` advance the reference, `N` is an intron, `I S H` do not; at least 2 inside junctions).
  Pools, counted as DISTINCT READ NAMES per member and chain (alignment class by flag: secondary if flag & 256, else supplementary if flag & 2048, else primary): **P1** = the alignment is primary (flag & 2308 == 0); **P2** = P1 plus secondary alignments (flag 256) whose read's primary alignment was found among the alignments fetched for the family's member spans
  and whose AS >= 0.98 x the best AS of that read (max AS over its non-supplementary alignments found there; a missing AS counts 0); **P2u** = secondary alignments carrying the chain whose read's primary alignment was NOT found in the family's spans (reported, never counted in P2, never gating); **P3** = any alignment, supplementary included.
  A member's support n = max over its valid chains of |P2 reads of that chain| (the support chain: ties by more introns, then lexicographic). **PASS iff n >= 3 for EVERY member.** A member without a valid chain has n = 0.
  Paralog-tie fraction of a member = (P2 reads of its support chain whose alignment at the member is NOT a primary with MAPQ > 0, i.e. secondary or MAPQ 0) / n; also printed: the share of those whose primary alignment lies on a different member of the family. Family tie fraction = pooled over members (sum tied / sum n).
  Descriptive, never gating: n_P2 with 5 bp junction tolerance (same number of inside junctions, every junction end within 5 bp of the chain).

## Overall pass, order of the counted steps, ranking

PASS-ALL = C1 and C2 and C3 and C4 and C5. Filter counts are reported cumulatively in the order C1 -> C2 -> C3 -> C4 -> C5 from N0 = 377 (the whole-projection step C1a is also printed as an intermediate line), and as independent per-criterion pass counts among the families on which the criterion is evaluable
(C1: all 377; C2, C3, C5: every member has a gene, i.e. the whole families; C4: every member has an RT).
Shortlist (class K) = PASS-ALL families ranked by (1) minimum over members of the P2 support n, descending; (2) minimum over members of the P1 count of the support chain, descending; (3) fewer members; (4) family number ascending. |K| is printed per size (2 / 3 / 4).
Near-miss list (descriptive): `near_miss.tsv` = families that fail exactly one criterion and pass the other four (the human code's rule); `near_miss_with_na.tsv` = exactly one failing criterion and no other failing one, with at least one criterion not evaluable (e.g. C1 fails for a missing RT, so C4 is NA).
A rule above is never edited after the first count; a code bug fix that changes a value is logged in `amendments_gorilla.md` with the reason and the before/after. This file is not edited after its hash is taken.

## Descriptive supplements (computed after the freeze, never gates or ranks)

- For every PASS-ALL family: the good-secondary cross-link (share of primary alignments in member A's span whose read has a secondary alignment in member B's span with AS >= 0.98 x the read's best AS; also the any-secondary and the >= 2-junction counts) exactly as the human `cross_links.py`, and the family tie fraction of C5.
- Truth-side columns copied, not recomputed: the number of protein-coding gorilla genes by name beyond the listed size (unlisted relatives) and LOC-named copies, per family.
- Genome-wide sanity: the gene-record coordinates in the truth table against the GFF parse.

## Independent re-computation

For every whole family in which every member has an RT (all of them if time allows, in chunks; in any case all PASS-ALL families, then all near-misses, then the lowest family numbers, a deterministic order; at least 10 families), C1 (RT, exon count, motifs, chains), C2 (shared bp), C4 (identity, coverage) and C5 (P1, P2, P2u, P3, tie) are recomputed by separate code that uses only command-line tools
(awk over the raw GFF, samtools faidx, bedtools merge / intersect, minimap2, samtools view text parsing), no pysam and no shared parser. Every mismatch is listed in `verification_gorilla.log`; a mismatch that comes from a code bug is handled by the amendment rule above.

## Outputs (all under SCRATCH)

`criteria_gorilla.md` (this file) + `criteria_gorilla.sha256` (hash and UTC time, taken before the first count), `inputs.tsv` (one row per family of the 377: every value and pass/fail per criterion, rank), `members.tsv` (one row per member of the 864), `pairs.tsv` (one row per member pair of the whole families),
`unmapped.tsv` (members that could not be projected or mapped, with the reason), `partial_families.tsv`, `shortlist.tsv`, `near_miss.tsv`, `near_miss_with_na.tsv`, `counts.json`, `cross_links.tsv`, `provenance.txt`, `amendments_gorilla.md`, `verification_gorilla.log`, the scripts, and the minimap2 FASTA / PAF files.

## Declared limitations of the design (not gates)

(1) the family size is the number of listed human genes, not the number of gorilla members: lineage-specific duplicates and LOC-named copies are invisible to the gate (the truth reports 72 of 377 families with more protein-coding gorilla genes by name than the listed size); (2) the truth is a naming layer (exact symbol, NCBI nomenclature) over a human gene tree: it over-represents older, dispersed pairs, and recent tandem duplicates that NCBI leaves LOC-named are absent (symbol-whole is 54% of cross-contig and 18% of same-contig families, quoted from the truth summary);
(3) reads (OR6737) and the assembly / annotation (KB3781) come from different individuals, so exact-chain matching can lose reads to individual variants at splice sites and copy-number differences are not modelled; chrX / chrY members are hemizygous in a male; (4) P2 resolves a secondary alignment's primary only inside the family's own spans, so P2 can be a lower bound (P2u shows how much);
(5) exact-chain matching ignores reads with 5' truncation or boundary jitter (the tolerance column shows how much); (6) the RefSeq annotation is Gnomon-heavy: many chains are model chains (XM_) that no read confirms; (7) one testis library, gorilla only; (8) the selection builds a positive-heavy set, so an O1 score on it is a mechanism check, not a recall estimate; the shortlist is not split into dev and held-out here;
(9) no genome-wide homology search: unannotated or non-coding paralogs elsewhere are invisible; (10) gorilla and human families share family-number labels only.
