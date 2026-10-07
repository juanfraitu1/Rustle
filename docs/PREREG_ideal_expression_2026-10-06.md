# PREREG: if every copy of NPIP and of TBC1D3 were ideally expressed, would the current default pipeline find them? (2026-10-06)

DRAFT v0, not yet committed. Written before any product of this protocol exists: no ideal read set of S1 has been simulated or mapped, and no S2 or S1 read set has been assembled with the HEAD binaries.
Everything here is DEV and human; the two families are never pooled. Nothing is tuned.

## Question

User, 2026-10-06: "like we tried earlier, if we had ideal cases of NPIP and TBC1D3 all expressed, would we be able to find them?" The earlier ideal experiments (`bench/NPIP_IDEAL_EXPRESSION.md` 2026-09-17, `bench/NPIP_SIM_CEILING.md` 2026-09-18,
`docs/IDEAL_CHROMOSOME_SIM_2026-09-21.md`) ran before f1v2, `--min-cov-shorter 0.70`, GOOD secondary seeding, the families copy table and the CAT default annotation, and the 09-17 reads were 30 identical copies per transcript (no end jitter: they collapse under the
assembler's coordinate dedup). They measured an algorithmic ceiling: with full-length reads and no readthrough NPIP had 25/26 complete chains (assembler perfect on what the reads contain), 27/27 copies had a node but only 11/27 a full-length one and 14/46 member nodes swallowed
a neighbouring locus; TBC1D3 was near-perfect (19/19 nodes, 17 full-length). This protocol repeats the question with the current default pipeline.

**Ideal** = every annotated transcript of every gene in the neighbourhood is expressed at equal depth, full length (ends jittered 0-30 bp), HiFi error only (substitution .001, indel .0003), no readthrough, no 5' truncation. It is CIRCULAR by construction (the reads come from the annotation that scores them): it
measures whether the default pipeline can find the copies when expression, truncation, readthrough and depth do not limit, not whether it recovers the annotation from real data.

## Substrates

- **S1 (new, CAT/Liftoff v2.0, the default annotation).** For each family, windows of +-500 kb around the territory (`terr_lo0`-`terr_hi`) of every CAT copy in `copy_recovery_tools_cat/ann/copies.hsa.tsv` (NPIP: 25 copies on chr16; TBC1D3: 16 copies on chr17), merged.
  Every CAT transcript of every gene (any biotype) whose span overlaps a window is simulated: 10 reads per transcript (a transcript with spliced length < 120 bp is skipped), read model of `bench/sim.py` (`simulate_reads`, err .001, indel .0003, 0-30 bp end jitter, seed 20261006),
  named `transcript|gene|fl|k`. The two families are simulated and mapped in separate read sets (no pooling). Mapping: the shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) against the whole CHM13 v2.0 splice index, in read-disjoint parts.
  **Canonicalization (the 2026-09-18 recipe):** an intron of a simulated transcript is kept only if it is >= 50 bp and canonical (GT-AG, GC-AG or AT-AC on the transcript's strand, checked on the genome); every other gap is merged into its flanking exons (the intronic sequence stays in the molecule).
  The canonical chain of a transcript is its remaining introns; transcripts with the same canonical chain on the same strand are one distinct chain of their gene.
- **S2 (existing, RefSeq, whole chr16; NPIP only).** `/mnt/linuxdisk/tmp/idealsim/ideal.bam` (4,411 RefSeq transcripts of 1,443 genes, 10 reads per transcript, jittered, err .001; `rt.bam` = the same plus readthrough molecules at the measured 7.51%, 3,014 of 43,700 reads) is assembled with the HEAD default pipeline,
  nothing re-simulated. The reads carry the annotation's own non-canonical junctions, which the default (`--assembly-junctions strict`) cannot assemble: a transcript is **recoverable** iff all its introns are >= 50 bp and canonical, and chain endpoints below use recoverable transcripts only
  (the number of copies with no recoverable transcript is reported: an annotation ceiling). RefSeq NPIP copies: `copy_recovery_tools/ann/copies.hsa.tsv`. TBC1D3 has no S2.

## Pipeline

`tools/rustle_pipeline.sh assemble` then `families` with the HEAD release binaries, no override, `--bam` the mapped reads, `--fasta` CHM13 v2.0 (secondary seeding from the best-AS table of that BAM, f1v2, shipped polish, `--min-cov-shorter 0.70`, most-reads representative). Products used: `PREFIX.gtf`
(all assembled transcripts), `PREFIX.families.gtf` (the loci the families stage read), `PREFIX.fam.clusters.tsv`, `PREFIX.fam.loci.gff3` (one representative transcript per locus), `PREFIX.fam.copies.tsv`.

## Targets and truth (fixed before any product is read)

- **Target copies:** S1: the 25 NPIP and the 16 TBC1D3 CAT copies; S2: the RefSeq NPIP copies of `copies.hsa.tsv`. A copy's truth = its gene's canonical chains (S1) / recoverable chains (S2) and the union of their exons.
- **Entangled copy:** another simulated gene on the same strand shares >= 100 exonic bp with the copy. Two annotated genes that share exons cannot get one locus each from any assembler; entangled copies are a truth ceiling, listed and reported beside, never removed from the denominator of the rules below.
- **Locus** = a gene_id group of `PREFIX.families.gtf`; its **representative** is the transcript of `PREFIX.fam.loci.gff3`. **Holder** of a copy = the same-strand locus whose representative exons overlap the copy's exon union by the most bp (ties: the larger locus).
  Coordinates are 0-based half-open exon blocks; chains compare intron (donor, acceptor) pairs exactly, ends are not compared.

## Endpoints (fixed denominators: all target copies of the family and substrate)

- **E0 observability** (descriptive): per copy, the share of its canonical junctions carried by >= 3 primary reads (-F 2308) of any MAPQ, and by >= 3 of MAPQ >= 1; reads of its transcripts that map back overlapping their source gene (primary).
- **E1 chains** (descriptive): per copy, distinct canonical chains simulated, and how many appear exactly as a transcript chain in `PREFIX.gtf`.
- **E2 own locus:** the holder exists, is the holder of no other target copy of the family, and its representative exons lie >= 50% inside the copy's exon union (not swallowing a neighbour: purity). Also reported: copies covered by >= 2 loci, copies whose holder is on the wrong strand.
- **E3 exact representative:** the holder's representative intron chain equals one of the copy's chains (S1: canonical; S2: recoverable). Relaxed E3': the representative's chain inside the copy's span is a contiguous sub-chain of >= min(2, n) introns of one of the copy's chains.
- **E4 in the family:** the holder is a member of the family cluster that holds the most target copies' holders. Also reported: the cluster's member loci that overlap no target copy (foreign members), and `family_score` (pooled bipartite sens / prec / F, pairwise) on the universe of simulated genes against Compara, U2 (chr16) and Soto.
- **FOUND** = E2 and E3 and E4. **FOUND-relaxed** = E2 and E3' and E4. **ENTANGLEMENT-FREE FOUND** = FOUND over the non-entangled copies (beside).

## Rules (per family and substrate)

- **YES, found:** FOUND >= 90% of the copies (NPIP 23 of 25 on S1; TBC1D3 15 of 16; S2 NPIP >= 90% of its copies).
- **PARTLY:** E2 and E4 each >= 90% but FOUND < 90% (the loci exist and are in the family, the representative structure falls short).
- **NO:** E2 < 90% or E4 < 90%.
- No rescue: no arm is re-run with another setting after a number is seen.

## Predicted before looking

- P1: TBC1D3 on S1 is YES (09-17: 19/19 nodes, 17 full-length).
- P2: NPIP on S1 is PARTLY or NO: loci exist for nearly all copies but some swallow neighbours or split (09-17: 14/46 member nodes swallowed a neighbour), and FOUND < 90%.
- P3: NPIP on S2 (whole chromosome, RefSeq) falls within 3 copies of NPIP on S1 in FOUND.
- P4: readthrough (S2 `rt.bam`) lowers NPIP FOUND by >= 3 copies against S2 ideal.

## Checks (a failure makes the arm INVALID)

G1 every simulated S1 chain is canonical and >= 50 bp (asserted from the genome); G2 >= 99% of the reads of the target copies' transcripts have a primary alignment overlapping their source gene; G3 the driver stages exit 0 and `families.gtf` is newer than `gtf`;
G4 the scorer is run twice with different `PYTHONHASHSEED` and gives identical tables; G5 an independent scorer written from this text alone reproduces E2, E3, E4 and FOUND per copy.

## Not claimed

Anything about real reads, the annotation's own correctness, O2 or O3, gorilla, or any family other than these two. A YES means the default pipeline has no algorithmic obstacle at these loci when the data are ideal; it does not predict recovery from a real library
(5' truncation cost 6 of 26 NPIP copies in 09-18, and real expression differs).
