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
