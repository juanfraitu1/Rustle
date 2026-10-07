# PREREG: if every copy of NPIP and of TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

Version 1, after an independent four-lens review of the first draft (metric traps, feasibility, truth tables, prior art; workflow wf_44161cf9-f29). Registered before any S1 read is simulated or mapped.
Human only, DEV, the two families never pooled, nothing tuned. **Disclosures** (what was read before this registration): (a) the first draft's smoke test of the scorer on the REAL default products (real reads, not ideal): NPIP E2/E3/E4/FOUND 11/3/20/3 of 25, TBC1D3 12/9/11/9 of 16;
(b) the reviewers computed truth-table statistics from the annotation (entanglement, chainless copies, canonicalization effects) and E0-like read statistics on the legacy S2 BAMs; (c) the HEAD default pipeline was run on the two legacy S2 BAMs (24 s and 71 s) and NOT scored; (d) the 2026-09-21 products on those BAMs were read when the first draft was written.
The first draft's rules were changed because the reviewers showed that its YES bars could not be met by construction and that its validity gate (G2) fails by design on this substrate; both are recorded below as the reason for the strata and for the new G2.

## Question

User, 2026-10-06: "like we tried earlier, if we had ideal cases of NPIP and TBC1D3 all expressed, would we be able to find them?" **Ideal** = every annotated transcript of every gene in the neighbourhood is expressed at equal depth, full length (ends jittered 0-30 bp), HiFi error only
(substitution .001, indel .00033), no injected readthrough, no 5' truncation. It is CIRCULAR by construction (the reads come from the annotation that scores them): it asks whether the default pipeline can find the copies when expression, truncation, readthrough and depth do not limit.

Earlier evidence, all older than the current default (f1v2, `--min-cov-shorter 0.70`, GOOD seeding, the families copy table, CAT annotation) and none a human run of the current default on ideal NPIP reads: 09-17 human nodes from 30 identical reads per transcript (NPIP 27/27 nodes, 11/27 full-length, 13/27 swallowing a neighbour; TBC1D3 19/19 nodes, 17/19 full-length;
node level only, never chain level for TBC1D3); 09-18 NPIP 25/26 complete chains (one canonical transcript per copy, no neighbours); 09-21 whole-chr16 (see "Legacy substrate"); gorilla 09-28 NPIP 25/25 complete and in the family with the copies alone, 24/25 in the family with neighbours, 22-24/25 under the current default at f = 0.

## Substrate S1 (new; CAT/Liftoff v2.0, the default annotation)

- **Windows:** +-500 kb around the territory (`terr_lo0`-`terr_hi`) of every copy in `copy_recovery_tools_cat/ann/copies.hsa.tsv`, merged: NPIP 25 copies on chr16 (7 windows, 12.19 Mb), TBC1D3 16 copies on chr17 (9 windows, 9.53 Mb). RefSeq NPIPB3 has no CAT/Liftoff gene, so S1 has 25 NPIP copies. 8 of the 41 target records are Liftoff, not CAT;
  everything is keyed by copy id (`cid`) and gene ID, never by name.
- **Transcripts:** every transcript of every gene overlapping a window (expected: NPIP 547 genes, 2,147 transcripts; TBC1D3 519 genes, 1,725). **Canonicalization (the 2026-09-18 recipe):** an intron is kept only if >= 50 bp and canonical (GT-AG, GC-AG, AT-AC on the transcript's strand, read from the genome); every other gap is merged into its flanking exons.
  A transcript is **simulated** iff its spliced length is >= 120 bp, its molecule is <= 30 kb and no kept intron exceeds 200 kb (minimap2's default maximum intron length, so such a molecule cannot align spliced); skipped transcripts are listed with the reason (expected: 77 and 55 under 120 bp, one 121 kb molecule, the transcripts with the 1 Mb intron of lncRNA AC009093.11 and one 218 kb intron).
  Annotated fusion and readthrough gene models are simulated like any gene (the arm has no INJECTED readthrough).
- **Reads:** 10 per simulated transcript, the model of `bench/sim.py` (`simulate_reads`: err .001, indel err/3, no truncation; 0-30 bp end jitter, skipped when the body would fall under 100 bp), each read drawn from `stable_seed(seed, transcript id, k)`. Two independent replicates, seeds 20261006 and 20261007, each simulated, mapped and scored.
- **Mapping:** the shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) against the whole CHM13 v2.0 splice index, in read-disjoint parts.
- **Families are simulated and run separately** (one read set per family and replicate: 4 mapped sets).

## Legacy substrate S2 (`/mnt/linuxdisk/tmp/idealsim/{ideal,rt}.bam`): NOT USED

Found by the review and confirmed here on 300 of 300 testable transcripts: `sim.py`'s `load_transcripts` keeps an exon only if its transcript row came earlier in the file, so every transcript's leftmost annotated exon (6,950 of 80,664 chr16 exon rows) is missing from the 2026-09-21 reads, which are therefore not full length;
the reads were mapped to `chr16.fa` only (one `@SQ` line, checked); and the `rt` arm replaces 7% of the reads and, by the reviewer's check of the reads, puts the downstream gene first on the minus strand (1,070 of its 1,109 minus-strand readthrough reads, 36.8% of all readthrough reads). The 2026-09-21 ideal-chromosome ceilings (`docs/IDEAL_CHROMOSOME_SIM_2026-09-21.md`) describe those reads. Scoring them would need its own registration.

## Pipeline and provenance

`tools/rustle_pipeline.sh assemble` then `families`, HEAD release binaries, no override (secondary seeding from the best-AS table of the same BAM, f1v2, shipped polish, strict junctions, `--min-cov-shorter 0.70`, most-reads representative). Products: `PREFIX.gtf`, `PREFIX.families.gtf`, `PREFIX.fam.loci.gff3` (one representative per locus), `PREFIX.fam.loci.tsv` (folded span to representative span),
`PREFIX.fam.copies.tsv`, `PREFIX.fam.clusters.tsv`. **G0:** sha1 of `copy_assign` (87824d91), `mcl_families` (a6308244), `as_table`, `family_score` (542923fd), the BAMs and the scripts is logged; the binaries must equal those recorded in `docs/DEFAULT_RESCORE_NPIP_2026-10-06.md`.

## Targets, strata and truth (fixed from inputs before any read is mapped; listed in `targets.tsv`)

Coordinates are 0-based half-open; an intron is (end of the upstream exon, start of the downstream exon); chains are tuples of introns in genomic order; one convention is asserted on a known intron (NPIPA1's first intron, 14938572-14942865).

- **E, entangled:** another simulated gene on the same strand shares >= 100 exonic bp (canonical exon unions) with the copy. The shared bp and the share of the copy's exon union are recorded; copies sharing < 50% are reported separately. Two annotated genes that share exons cannot each get a locus from any assembler.
- **C, chainless (MONO):** the copy has no multi-exon canonical chain (NPIPA3 and TBC1D3P7: every raw intron is non-canonical, so the molecules are unspliced gene bodies). Unspliced loci default to the + strand, so a MONO copy's holder is strand-agnostic.
- **X, readthrough image:** the copy table has `readthrough = 1` (NPIP h03, h05, h08: images of fusion records).
- **R_in = copies in none of E, C, X.** Expected before E0: NPIP <= 14 of 25, TBC1D3 <= 13 of 16. **R** (per replicate) = R_in copies with an OBSERVABLE chain (E0 below); R is computed from the truth tables and the BAM alone, before assembly, and its sha1 is logged.
- **Named-family genes:** simulated non-target genes whose CAT/Liftoff `gene_name` contains NPIP or TBC1D3 (e.g. NPIPP1, PKD1P6-NPIPP1, NPIPB7, PDXDC2P-NPIPB14P, TBC1D3J, TBC1D3P1-DHX40P1): expected non-target family members, never "foreign".
- **Single-copy genes** (instrument controls): simulated genes none of whose reads has a secondary alignment with AS >= 0.9 x its best AS, taken from the BAM.

## Endpoints (denominators fixed from inputs)

- **E0 observability (BAM alone).** Pools: P1 = primary alignments (-F 2308); P2 = P1 + secondary alignments with AS >= 0.98 x the read's best AS (the pipeline's seeding pool); P3 = any alignment. A chain is OBSERVABLE iff >= 3 reads of P2 carry exactly that whole chain (the read's junctions inside the copy's span equal the chain, N >= 50 bp).
  Printed per copy and pool: reads, share whose primary overlaps the source copy (no bar), observable chains. A copy with no observable chain is ALIGNER-LIMITED.
- **Holder.** Records are resolved through `fam.loci.tsv` (folded loci count as their representative; the number of folded holders is reported). The holder of a copy is the same-strand locus maximizing J = overlap(rep exons, copy exon union) / union of the two (bp); ties: fewer representative exonic bp, then lower gene name.
  Clusters join to loci by (chrom, gene-row 1-based start, end) as `bench/default_rescore/nodes.py` does; two loci on one span or an unmatched cluster row make the arm INVALID.
- **E1 chains:** distinct canonical multi-exon chains of the copy that appear exactly as a transcript chain in `PREFIX.gtf`, over all chains and over observable chains (descriptive).
- **E2 own locus:** a holder exists, holds no other target copy of the family, and both purities are >= 0.5: representative purity |rep exons in the copy| / |rep exons| and locus purity |exon union of the locus's transcripts in the copy| / |that union|. Also printed: same-strand loci overlapping >= 100 exonic bp (>= 2 = split), and whether the top-overlap locus is on the other strand.
- **E3 exact representative:** the holder's representative chain equals one of the copy's chains. **E3':** the representative's introns inside the copy span are a contiguous sub-chain of ONE chain with >= max(2, ceil(n/2)) introns, n = introns of that chain. **E3c** = chains recovered exactly among the holder locus's transcripts / chains simulated; **E3p** = holder-locus transcripts with a chain that is not a chain of the copy / holder-locus transcripts (over-emission).
  MONO copies: E3 = the representative covers >= 90% of the copy and has <= 110% of its bp (reported beside, outside R).
- **E4 in the family:** K* = the cluster of the holder of the lowest-cid R copy that has a holder in a cluster; the holder is in K*. **CP** = (member loci of K* overlapping a target copy or a named-family gene) / |K*|. Per family: pairwise tp / truth / predicted pairs, then bipartite sens / prec / F with `family_score` on truth rows restricted to the simulated windows, predicted members outside the truth counted as false, with and without size-2 clusters; the pooled line is printed labelled not-a-family-result.
- **IDEAL-FOUND = E2 and E3 and E4** (this is NOT the registered strict FOUND of Amendment E, which needs the cap signal and cannot run on simulated reads). **IDEAL-FOUND-relaxed** uses E3'. Printed beside, on the same copies: the `nodes.py` own-node flag and `bench/copy_support.py` Amendment A `ann_found` / `chain_found` computed on the ideal BAM, next to the registered real-read DEF values
  (own node 24/25, strict FOUND 8/25, locus level 21/25, U2 F .645, CF153 13/19 F .812); only own node and `ann_found` are like-for-like.

## Rules (per family and replicate, on R; n = |R|, N = all target copies of the family)

- **CEILING-LIMITED** if n < 0.5 N: no YES or NO is issued.
- **YES:** IDEAL-FOUND >= ceil(0.9 n) and CP >= 0.5.
- **PARTLY:** E2 and E4 each >= ceil(0.9 n) but not YES.
- **NO:** otherwise.
- The family's verdict is the LOWER of the two replicates (NO < PARTLY < YES); if the replicates differ it is reported UNSTABLE and no YES is claimed. The copies whose flags differ between replicates are listed.
- Always printed beside, never instead: IDEAL-FOUND over ALL copies as x/N, and one line per excluded stratum (E, C, X, ALIGNER-LIMITED). An arm that fails G1, G3, G4 or G6 is INVALID; an INVALID arm is amended in writing before any pipeline product of that arm exists, afterwards only a new arm may be registered. No rescue.

## Controls

- **G1** every simulated chain is canonical and >= 50 bp, re-read from the genome by a separate script (`verify_g1.py`) before mapping. **G2 (instrument)** >= 99% of the reads of single-copy genes have a primary alignment overlapping their source gene; the same fraction for the target copies is E0, printed, ungated (near-identical paralogs coin-toss primaries).
  **G3** the BAM holds exactly one primary-or-unmapped record per FASTQ read; the driver stages exit 0 and `families.gtf` is newer than `gtf`. **G4** the scorer is run twice under different `PYTHONHASHSEED` values and gives identical tables.
- **G5** an independent scorer written from this text alone (by someone who has not seen `score.py`) reproduces E2, E3, E4 and IDEAL-FOUND per copy.
- **G6 positive control:** the canonicalized annotation as loci (each gene a locus, equal read counts) through `mcl_families --from-gtf` with the shipped flags and the same scorer: E2 and E3 must be 100% on R, and its E4 is the family ceiling E4*; YES needs E4 >= min(0.9, E4*) on R.
  **G7 negative control:** >= 95% of single-copy simulated genes get exactly one locus (same strand, >= 100 exonic bp) and none lies in K*.

## Held-out contrast (no bars, labelled HELD-OUT: the pipeline was developed on NPIP and TBC1D3)

Two other multi-copy families inside the same windows, chosen now by CAT/Liftoff `gene_name`: KRTAP genes (`^KRTAP[0-9]+-[0-9]+$`) in the TBC1D3 windows and SMG1P genes (`^SMG1P[0-9]+$`) in the NPIP windows. Scored with E2-E4 on their non-entangled copies that have a canonical multi-exon chain; membership is not verified before the run.

## Predicted before looking

- **P1** TBC1D3 is YES on R in both replicates (09-17: 19/19 nodes, 17/19 full-length).
- **P2** NPIP is not YES on R (09-17: 13 of 27 nodes swallowed a neighbour) and IDEAL-FOUND over all 25 copies is <= 21.
- **P3** the two replicates give the same verdict in both families.

## Also reported with every headline

Per-copy tables (E0 three pools, E1, E2 with both purities and the holder, E3/E3'/E3c/E3p, E4, entangled-with gene, shared bp and share, stratum, and one miss reason from the fixed order annotation, aligner, locus formation, family membership, representative); the seed spread; per family CP with every non-target member by gene name and biotype;
the single-copy control; chain-level recall and precision; binary and BAM sha1. One sentence of scope: DEV families, human only, circular by construction, +-500 kb neighbourhood only (precision against homologous loci elsewhere is not measured), equal-depth full-length reads, annotated fusion records present, Liftoff records among the targets.

## Not claimed

Anything about real reads, the annotation's correctness, O2, O3, gorilla, or any family beyond these two and the two held-out ones. A YES means the default pipeline has no algorithmic obstacle at these loci when the data are ideal; it does not predict recovery from a real library (5' truncation cost 6 of 26 NPIP copies in 09-18, and real expression differs).

## Amendment 1 (2026-10-06, before any S1 read exists)

The shared-span rule of the holder paragraph is narrowed: two loci on one (chrom, start, end) span make an arm INVALID only when a cluster row or a fold row refers to that span (the join would be ambiguous); shared spans among loci that no cluster or fold row names are ignored.
The first scorer test on the REAL default chr17 products showed two unclustered loci on one span (`DN_chr17_87073_2` and `_3`), which the draft rule would have declared invalid although no cluster or fold row names them. Nothing else changes.

