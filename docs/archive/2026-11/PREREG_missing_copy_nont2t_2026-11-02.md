# PREREG MISSING COPIES ON A NON-T2T GENOME (register 1504, 2026-11-02, before any alignment or number)

Question (user 2026-11-02, after the advisor's "a T2T genome and a good annotation will not exist for every organism"): with the
same Iso-Seq reads aligned to an incomplete assembly, do our O3 stages find the copies the T2T assembly has and the incomplete one
lacks? This run is the BASELINE (today's tools, no tuning); the read signals it records are what later rules learn from.
Pairs (same species, reads already aligned to the T2T): human A119b, CHM13 (T2T) vs GRCh38 no-alt analysis set; gorilla GGO_mm,
mGorGor1 (T2T) vs Kamilah_GGO_v0 (GCF_008122165.1). Downloads in ~/annot/nont2t/.
Panel families: the famtruth MAINSTAYS and YAGS patterns (NPIP, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE, SPATA31; BPY2, CDY, DAZ,
HSFY, PRY, RBMY, TSPY, VCY, XKRY) on the T2T annotation's gene names (human CAT; gorilla mGorGor1 RefSeq, plus GGO_MAINSTAYS).
Reads: every read name with any record (primary or secondary) on a panel gene of the T2T BAM, re-aligned with the T2T BAM's own
minimap2 settings (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) to the non-T2T genome.
Truth (sequence only, no tool of ours): each panel T2T copy's spliced exon sequence (longest isoform) is aligned to the non-T2T
genome (minimap2 -x splice:hq, all hits). A copy is PRESENT when some hit has identity >= 99% over >= 95% of the copy and that
locus (hit span, same strand) is not already claimed by a better-identity copy; else MISSING (the incomplete genome lacks or
collapsed it). EXPRESSED = >= 3 primary reads overlapping the copy's exons in the T2T BAM. Truth set = expressed missing copies.
Arm: on the non-T2T BAM, `rustle_pipeline.sh assemble`, `families`, `candidates`, then `flag` (no --confirm, no --foreign: no T2T
leak). Calls: flag rows with verdict `reference_absent_candidate`; candidates rows with `flagged` true.
Matching (scoring only): each call's consensus / candidate sequence is aligned to the T2T panel copies; a call HITS a missing copy
when its best T2T copy is that copy at identity >= 99%.
Report per pair and per tool (flag, candidates, union): recall = truth copies hit / truth copies; precision = calls hitting any
missing copy / calls; calls hitting a PRESENT copy (false alarms on copies the genome has) listed separately.
Decision: none on defaults (baseline). Pre-declared use: if union recall >= 0.5 and precision >= 0.5 on BOTH pairs the tools are
reported as usable on incomplete genomes as they are; otherwise rule learning follows with HUMAN as development and GORILLA held out
(rules fixed on human, scored once on gorilla).

## Amendment 1 (2026-11-02, before any run)
A copy's sequence = the spliced exon union of its gene (famtruth.load_cat for human, load_ggo for gorilla: the loaders that carry
the gene names), not the longest isoform (the CAT slim GFF3 has no gene names). Same for both pairs.

## RESULT r1504 (2026-11-02): FAIL — neither tool finds the missing copies
Job 4417876; `/anvme/workspace/v138da13-rustle/nont2t/{human,gorilla}` (`score.py`, calls at >= 99% identity to a T2T panel copy).
Flag calls = verdicts `fired*` (`x.flag.missing_copy.tsv`); candidates = `flagged == 1`.
| pair | arm | calls | missing expressed | found | recall | precision | on present copies |
|---|---|---|---|---|---|---|---|
| A119b -> GRCh38 | flag | 63 | 32 | 1 (NPIPB4) | 0.031 | 0.016 | 5 |
| A119b -> GRCh38 | candidates | 100 | 32 | 1 (NPIPB4) | 0.031 | 0.010 | 4 |
| A119b -> GRCh38 | union | 163 | 32 | 1 | 0.031 | 0.012 | 9 |
| GGO -> Kamilah | flag | 21 | 37 | 0 | 0 | 0 | 1 |
| GGO -> Kamilah | candidates | 44 | 37 | 0 | 0 | 0 | 2 |
| GGO -> Kamilah | union | 65 | 37 | 0 | 0 | 0 | 3 |
Not a threshold effect: at >= 95% identity human stays 1/32, gorilla reaches 5/37 (precision 5/65). Most calls resemble no
panel copy: human 55/163 have no hit at all and 71 best-hit a copy below 95% (median query coverage of the best hit 0.49);
gorilla 42/65 below 95%. The calls are fragments and mixtures of present loci, not the T2T copies the genome lacks (TSPY, GOLGA6L,
NPIPB, DAZ, LOC129* families). The pre-declared else-branch (rule learning, human dev / gorilla held out) has nothing to learn from
at recall <= 0.03: not run. Verdict: the flag and candidates tools are NOT usable as missing-copy detectors on a non-T2T genome
as they stand; a detector would have to start from the reads that fail to place (the `~xchrom~`/unplaced net) rather than from
mixtures at present loci. Timings (72 cores): human realign 660 s, assemble 68 s, families 4 s, candidates 1,456 s
(candidate_copies 834 s: nets 357, clusters 382, genome 25, structure 69), flag 403 s; gorilla 128 / 21 / 6 / 358 / 44 s.
