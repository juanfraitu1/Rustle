#!/bin/bash
# rustle_pipeline.sh — the whole pipeline, one command per stage, shipped defaults (2026-09-23).
#
#   assemble  reads -> loci -> isoform GTF                              copy_assign --assemble-only --genome-wide
#   families  gene families from the de novo loci (all-vs-all -> MCL)   mcl_families --from-gtf --emit-units
#             = THE default de novo family definition (user decision 2026-09-25): one copy per member locus =
#             its representative transcript (PREFIX.fam.copies.tsv/.fa, the copy table copy assignment consumes)
#   candidates  reference-absent copies from each family's own reads    candidate_copies + utilities candidate-augment + minimap2
#             OPT-IN (ruling R14, 2026-10-02: its pre-registered acceptance, Amendment 12, FAILED —
#             docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; 2026-10-03: the re-run, Amendment 13, passed, but the no-deletion
#             control, Amendment 14, failed, so the default flip was reverted — docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md):
#             naming the stage runs it, `all` runs it only with --candidates.
#             (read net -> read clusters -> consensus -> flag/link/merge -> one exon-union contig per flagged
#             candidate, PREFIX.cand.*), then the augmentation (PREFIX.aug.fa = genome + the contigs, PREFIX.aug.copies.*
#             = the copy table + one row per candidate) and the patch realignment of the candidate families' reads to
#             PREFIX.aug.fa (PREFIX.aug.bam): a read whose best alignment is on a candidate contig is PLACED there by the
#             realignment (ruling R13). --delta D (0.00958), --cand-max-reads N (1000). Cost: PREFIX.aug.fa is a copy of
#             the whole genome, and the realignment indexes it on every run (about 400 s and 19 GB on a human genome).
#             The stage's own cost on a full BAM (spec §9b, R23): on the 23-GB gorilla fibroblast BAM one batch of 50 families
#             did not finish in a 10-minute call (nets phase 397-469 s: the whole-BAM sweep ~150 s, the attribution alignment
#             against the batch's 2.6 GB of net reads 213-250 s; peak RSS 10.7 GB at the stop). The driver runs EVERY family in
#             one call (no batching; the binary's `--families` batches): the attribution targets are then all families' nets, so
#             its peak memory exceeds that batch's 10.7 GB — an OOM risk on this machine.
#   catalog   LEGACY copy catalog (gw_family_catalog; kept, not the default definition)
#   assign    per-read copy assignment on the catalog (assign/abstain)  copy_assign --families
#   all --to STAGE / --from STAGE   stop after / start at a stage (escape hatch; `--to families` answers "what do the loci
#             look like before reads are assigned": PREFIX.fam.copies.tsv + PREFIX.fam.loci.gff3). 2026-10-30.
#   merged    THE genome-wide entry point: assemble + families + assign on the FAMILIES copy table, one command,
#             byte-identical products, with the genome-wide levers ON by default (--skip-poa-diagnostic and
#             --region-threads, both cmp-proven byte-identical) and per-phase resume (a finished phase is never
#             redone; the families all-vs-all PAF replays from PREFIX.cache when the loci FASTA is unchanged).
#   assign    per-read copy assignment (assign/abstain) on the families'  copy_assign --families
#             copy table PREFIX.fam.copies.*; with --candidates and flagged candidates, two runs (the candidate families
#             on PREFIX.aug.*, each assigned as a family that includes its candidate copies; the others on the originals)
#             concatenated into PREFIX.assign.*. O2 assigns AS-tied molecules only (the gate is unchanged): a read the
#             realignment places uniquely on a candidate has no row (R13). Without --candidates the candidates stage's
#             products are not used. --legacy-catalog: on the legacy catalog PREFIX.cat.* instead (the stage as it ran
#             before 2026-10-02; refused together with --candidates)
#   flag      copies the reference does not contain, from RNA alone     missing_copy_flag --scan-only / --from-scan
#             (+ optional DNA confirmation against --confirm genomes; with --candidates, + the `o3_candidate` column
#             naming the flagged candidate whose nearest locus is the row's)
#   all       assemble, families, o2, flag (register 1903: `o2` in place of the legacy `assign`, in every plan);
#             with --candidates the candidates stage before o2, with --legacy-catalog the catalog there
#   sd        the genome's segmental duplications and duplicons            utilities sd-scan (vendored sedef-rs)
#             once per genome, cached by the FASTA's md5 and sd-scan's version key in RUSTLE_SD_CACHE_DIR (default
#             PREFIX.cache)/sd/<key>/genome.sd.*; PREFIX.sd.prefix names that prefix for `mcl_families --emit-hierarchy
#             --sd-cache`. Needs --fasta and --out only; `all` does not run it. An entry is a hit only if its DONE lists
#             every file with its size (an empty DONE, as this stage wrote before 2026-10-14, is a miss); --no-cache scans
#             again, into the same entry.
#             Gates in docs/archive/2026-11/LEDGER_O1.md:4038 (PREREG_sedef_rs_layers_2026-10-13.md) (register 1278-1280): sd-scan's SD calls failed V1/V1b
#             (recall about 0.40 against SEDEF) and its duplicons failed V2 (held-out family score 0.5153 against
#             DupMasker), so `--sd-cache` is an opt-in, unvalidated alternative to `--hierarchy-sd` /
#             `--hierarchy-duplicons`; its SD layer uses every `ok` pair (no identity or length floor). Cost with
#             --threads 4: about 18 min and 10 GiB peak per human genome, about 55 min and 19 GiB for gorilla. The stage
#             sets no memory cap (bench/sedef/run_genome.sh runs it under `ulimit -v`, 23 GiB).
# (thesis record: families/catalog = O1, candidates = O3, assign = O2, flag = O3; in order O1 families -> O3 candidates ->
#  O2 assign -> O3 flag)
# 2026-11-22 (user, register 1903): O2 is the family benchmark across tools on planted and real reads
# (bench/o2_families/, docs/O2_FAMILY_BENCHMARK.md). In `all` the stage `o2` REPLACES `assign`:
#   o2        what this run reports per gene family: copies, copies with a spliced transcript, transcripts, isoforms,
#             reads -> PREFIX.o2.families.tsv, PREFIX.o2.copies.tsv (tools/o2_family_report.py; seconds, no BAM pass);
#             with PREFIX.eq.loci.tsv (stage eqclasses) or PREFIX.multicopy.eq.loci.tsv present, also reads placed on
#             one copy, tie groups and the range of expressed copies. No read is assigned.
#             First, by DEFAULT since 2026-11-22 (user; RUSTLE_O2_EXTEND=0 turns it off): the assembly is extended to the
#             multi-copy families with the near-tied multimappers (= stage multicopy on the run's assembly, made with
#             --seed-pool good as `all` does -> PREFIX.multicopy.gtf), then the provenance test labels every locus by
#             the EM over the tied reads (bench/o2_families/em_filter.py; rows 1905-1908): o2_status = expressed /
#             shared / homology_only (under RUSTLE_O2_MIN_EM = 2 reads by EM) -> PREFIX.o2.gtf, PREFIX.o2.em.tsv.
#             Since 2026-11-23 (user, register 1921) the alignment scores are inside that EM: every alignment
#             within 10 ms of a read's best is a candidate with weight exp(-gap / RUSTLE_O2_EM_WEIGHT), default 1;
#             RUSTLE_O2_EM_WEIGHT=0 = the counts-only EM of before (scores decide exact ties only).
#             RUSTLE_O2_SEQUENCES=1 (opt-in, register 1926; needs --fasta): how many distinct sequences the tied reads
#             show among copies the reference cannot separate (IsoCon's test on genome alignments) ->
#             PREFIX.o2.sequences.tsv and the report column tied_sequences. Sequences, not copies: alleles count.
#             Its collected reads are cached in PREFIX.o2.sequences.cache (RUSTLE_O2_SEQUENCES_CACHE=0 turns that off).
#             An isoform with more than 1,000 distinct variant patterns (immunoglobulin, T-cell receptor) is not tested:
#             PREFIX.o2.sequences.skipped.tsv.
#             The pass over the BAM that finds, for every multimapping read, the loci it can sit at is made once
#             (tools/tie_table.py -> PREFIX.o2.ties.tsv.gz, readable: read, locus, ms gap, best ms) and read by the
#             labels and by the sequence test.
#             Nothing is dropped: a family member that is found stays available. RUSTLE_O2_DROP=1 drops the
#             homology_only loci instead. Costs a second assembly at a floor of one read, a families run on the
#             candidates and one BAM pass (gorilla OR6737: about 7 min). If it fails, the report is made on PREFIX.gtf
#             and the stage still succeeds. PREFIX.gtf itself is never changed. Needs python3 + pysam for the test.
#   o2bench   OPT-IN, never in `all`: the VALIDATION of O2, planted reads of known origin for every assembler on PATH
#             (bench/o2_families/plant.py + run_tools.sh) -> PREFIX.o2bench.tsv. With --gff ANNOT.gff3.gz: the named
#             families of the annotation (famtruth environment, human CHM13). Without --gff: this run's own families
#             (PREFIX.gtf + PREFIX.fam.copies.tsv, any genome; needs FASTA.fai). Scenarios RUSTLE_O2BENCH_SCENARIOS
#             (default E,D,S), reads per copy RUSTLE_O2BENCH_READS (default 60 named / 30 own families).
#   provenance OPT-IN, never in `all`: where each read went -> PREFIX.provenance.reads.tsv.gz (one row per primary alignment:
#             the transcript it supports and how: exact, after a junction replacement, as a shorter read, none and why;
#             short-exon rewrite, near tie, label, sequence), PREFIX.provenance.transcripts.tsv (the reads behind every
#             transcript), PREFIX.provenance.family_reads.tsv.gz. Reconstructed from the BAM and PREFIX.o2.gtf (or
#             PREFIX.gtf); what is inferred is marked in tools/read_provenance.py.
#   members   OPT-IN, never in `all` (register 1910): the relatives of every assembled locus through ALL of its reads'
#             alignments, not only the near ties (tools/member_discovery.py; RUSTLE_MEMBERS_MIN_READS 5,
#             RUSTLE_MEMBERS_MIN_RATIO 0.3) -> PREFIX.members.tsv (place, source locus, reads, score ratio, the assembled
#             locus there or `.`) and PREFIX.members.gtf (a model per place with no assembled locus, o2_status
#             "member_candidate"). Where relatives are, expressed or not. Needs a BAM mapped with the usual -p 0.1: the
#             slim mapping mode removes these records (row 1909). Needs python3 + pysam.
#   `assign` is the LEGACY O2 and is in NO plan of `all` any more (2026-11-22): not by default, not with --candidates,
#   not with --legacy-catalog, and there is no switch. The stage itself is kept and runs when called by name
#   (`rustle_pipeline.sh assign ...`, as tools/rustle_samples.sh does); the candidates and catalog stages still run in
#   `all` when asked for and leave their tables for such a call.
#
# usage: tools/rustle_pipeline.sh STAGE --bam B --fasta G --out PREFIX [--index G.splice.mmi] [--gff ANNOT.gff]
#        [--confirm NAME=X.mmi ...] [--foreign NAME=X.mmi ...] [--threads N] [--bin DIR]
#        [--seed-pool primary|good|all [--seed-as-ratio R]] [--as-table FILE] [--contig NAME[,NAME...]] [--no-seed-secondaries]
#        [--no-cache] [--inspect] [--piecewise [--max-pieces N] [--budget-s S] [--piece-records R] [--piece LABEL]]
#        [--candidates | --no-candidates] [--delta D] [--cand-max-reads N] [--legacy-catalog]
#        [--to STAGE] [--from STAGE] [--low-mem] [--index-gb N] [--max-window-records N] [--no-flag] [--no-resume] [--cache]
#        (stage `index` = the shared splice index alone, needs --fasta and --out; `flag_scan` / `flag_bg` are internal to `all`;
#         RUSTLE_EXTEND_TRANSFRAGS=1 adds --extend-transfrags to assemble, default off since register 1812)
#        [--guide ANNOT.gtf]   assemble / all also write PREFIX.guided.gtf (stage `guided` alone rewrites it from PREFIX.final.gtf;
#                              the de novo PREFIX.gtf and every later stage are unchanged; register 1622)
#        [--names ANNOT.gff3|gtf[.gz]]   assemble / all also write PREFIX.named.gtf (and PREFIX.final.named.gtf): the same
#                              transcripts with ref_gene_name / ref_gene_id / ref_transcript_id / ref_match attributes from the
#                              annotation (tools/gtf_names.py; exact intron chain, else the gene sharing most exon bases on the
#                              strand). Names only: independent of --guide, no product or stage reads the named files, ids and
#                              models are unchanged. Stage `names --names ANNOT --out PREFIX` alone rewrites them.
#   tools/rustle_pipeline.sh psi --bam B --fasta G --out PREFIX [--threads N]   OPT-IN, named only: PREFIX.psi.tsv = PSI of each
#     intron within its bunches (introns sharing a donor, or an acceptor; MntJULiP's unit) counted from the alignments, with the
#     range [psi_lo, psi_hi] that reads tied between placements allow, in eqclasses' two bands (tools/psi_bunches.py; needs
#     PREFIX.gtf and the score table). A report of one sample: no test between samples, intron retention not seen, nothing reads it.
#     After stage `exonmap` it also writes PREFIX.fam.psi.tsv: the psi rows keyed by exon columns, so one junction is read across
#     the copies of a family (tools/psi_family.py; `differs` = copies whose tie ranges do not overlap, a count, not a test).
#     With PREFIX.fam.loci.gated.paf it also adds rows for each copy that has the junction's two exons and reads at its anchor
#     but no row: its own junctions between those exons that the read floors dropped, else the junction carried over from the
#     most-read copy (`silent`: no read carries it), and the columns n_eligible, range_lo, range_hi, gap;
#     PREFIX.psi.introns.tsv = the raw counts of every intron, written for that purpose.
#   tools/rustle_pipeline.sh eqclasses --bam B --fasta G --out PREFIX [--names ANNOT] [--threads N]   OPT-IN, named only:
#     PREFIX.eq.classes.tsv (sets of loci minimap2's ms score does not separate, with their read count) and PREFIX.eq.loci.tsv (per
#     locus of the assembly: reads placed only there, the maximum with every tied read, a counts-only EM estimate between the
#     two, and its tie group). Quantification and reporting only: no read is assigned and no other product changes
#     (tools/eq_classes.py). Two bands side by side: exact score ties, and (columns _w) scores within RUSTLE_EQ_WINDOW of the
#     best (default 10 ms; or >= RUSTLE_EQ_RATIO x the best when that is set); neither is validated, so read where they disagree. Loci mostly
#     inside the exclude BED (RUSTLE_EXCLUDE_BED or FASTA.exclude.bed) are labelled repeat and join no group.
#   tools/rustle_pipeline.sh famwise --bam B --fasta G --out PREFIX [--names ANNOT] [--threads N]   OPT-IN, named only:
#     `eqclasses` plus PREFIX.famwise.*: transcripts of different loci that the same tied molecule carries are joined into
#     family transcripts (families.tsv: copy set, each molecule once); per transcript the molecules placed only there and
#     the maximum with every tied one (tx.tsv), and a status: certified, best_only, tied_only, near_only, no_band_record.
#     Bands on minimap2's ms (exact, and within RUSTLE_EQ_WINDOW); every record row also carries AS (reads.tsv, within 20).
#     Reporting only: no read is assigned, no transcript added or dropped, no other product changes (tools/famwise.py).
#     RUSTLE_FW_MIN_LINK (default 1) = tied molecules needed to join two transcripts; RUSTLE_EQ_RATIO is refused.
#   tools/rustle_pipeline.sh siblings --bam B --fasta G --out PREFIX   OPT-IN (register 1628): PREFIX.sibling.gtf = PREFIX.gtf +
#                              one-read chains a family sibling's model confirms; after assemble and families
#   tools/rustle_pipeline.sh multicopy --bam B --fasta G --out PREFIX   OPT-IN (registers 1858-1867, user 2026-11-21), the
#                              multi-copy mode; run after `assemble --seed-pool good` (near-tie secondaries seed the copies primaries
#                              miss). PREFIX.multicopy.gtf = PREFIX.gtf + spliced one-read transcripts at places PREFIX.gtf leaves
#                              empty, kept only when their locus is a copy of a gene family (tools/multicopy_rescue.py);
#                              PREFIX.multicopy.tsv lists them. Transcripts also carry counts that do not multiply over
#                              indistinguishable copies (reads_em, reads_tieadj, tie_group); PREFIX.multicopy.groups.tsv = the tie
#                              groups with the range of expressed copies (row 1871). No other product changes.
#   tools/rustle_pipeline.sh oneread --bam B --fasta G --out PREFIX [--threads N]   (registers 1766 / 1806 / 1808; also a last step of `all` since 2026-11-20,
#                              RUSTLE_ONEREAD_TIER=0 = off): side files, never read by another stage. PREFIX.oneread.gtf = one-read chains with every
#                              junction in another read and a clean reading frame (tools/oneread.py); PREFIX.oneread.ranked.gtf = every one-read chain
#                              the default lacks, scored (a ranking, not a filter). `--partners A.bam,B.bam` (other libraries of the species; =
#                              RUSTLE_ONEREAD_PARTNERS=A.bam:B.bam): PREFIX.oneread.promoted.gtf = ranked chains in the best 20% that a partner read
#                              carries, appended to PREFIX.final.gtf with oneread_promoted "1". RUSTLE_ONEREAD=1 also appends the strict tier,
#                              RUSTLE_ONEREAD=0 keeps PREFIX.final.gtf free of one-read models. A multi-sample run (--bam a,b) pools the libraries
#                              first, so a chain seen once in each already has two reads there. After assemble.
#                              RUSTLE_ONEREAD_GATE=0|1 (default 1 since row 1868, 2026-11-21; 0 = before): the candidates include the
#                              single-read chains the assembly gate's pooling guard hides (row 1864). The promoted tier is a RECALL tier
#                              (0.22-0.30 precision on three libraries, row 1866), below the default's two-read chains.
#   tools/rustle_pipeline.sh flairlike --bam B --fasta G --out PREFIX [--threads N]   OPT-IN (register 1640): PREFIX.flairlike.gtf =
#                              FLAIR collapse's chain set (second minimap2 pass; lower precision than PREFIX.gtf)
#   tools/rustle_pipeline.sh famreps --bam B --fasta G --out PREFIX   the second half of `families` alone (register 1662):
#       genes that have no locus of their own join PREFIX.fam.copies.* as extra representatives of their locus;
#       PREFIX.fam.base.copies.* keeps the table of the first run. RUSTLE_FAMREPS=0 turns it off inside `families`.
#   tools/rustle_pipeline.sh altends --bam B --fasta G --out PREFIX   OPT-IN: PREFIX.altends.gtf = PREFIX.gtf + alt_tes attributes
#                              (the commonest read 3' end when it is not the emitted one; polyA-signal and A-rich flags)
#   tools/rustle_pipeline.sh repeatloci --bam B --fasta G --out PREFIX   OPT-IN (registers 1698 / 1702): PREFIX.repeat_loci.tsv marks
#     assembly loci that are soft-masked repeats or rDNA (a table; nothing is dropped). Needs python3 + pysam + minimap2.
#   tools/rustle_pipeline.sh exonmap --bam B --fasta G --out PREFIX [--threads N]   OPT-IN (register 1595): PREFIX.fam.exonmap.columns.tsv =
#     the exon columns of each family, each core (every member), accessory (some) or private (one); PREFIX.fam.exonmap.gfa/ = a graph
#     view per family (segments tagged with the class, one path per member, Bandage colours). Reads the families output only.
#   tools/rustle_pipeline.sh relatives --bam B --fasta G --out PREFIX   OPT-IN (register 1676 / 1737): PREFIX.fam.relatives.tsv = locus pairs
#     the family edge rule refused that the row-1676 score calls related. Links, NOT family members; the families are unchanged.
#     Set RUSTLE_EXCLUDE_BED (the satellite BED assign uses): without it 3 of 4 links on A119b are satellite loci.
#     Needs a SEDEF BED (RUSTLE_SEDEF_BED or FASTA.sedef.bed). Measured only on protein-coding gene pairs: precision 0.830 on
#     two gorilla samples (its own bar 0.85 not passed), ~1/6 of the picks young duplicates; other locus pairs are unscored.
#   tools/rustle_pipeline.sh sdcore --bam B --fasta G --out PREFIX   (register 1755 / 1756; also the last step of `families` when the genome
#     has a SEDEF BED, RUSTLE_SDCORE=0 = off; default since 2026-11-15, user): PREFIX.fam.sdcore.tsv = pairs of copies of
#     one family that one SEDEF pair holds (the --sd-damp span rule) and aligns exon onto exon (sd_core = 1). A label; the families
#     are unchanged. Needs a SEDEF BED with CIGARs (RUSTLE_SEDEF_BED or FASTA.sedef.bed). HG002 held out: flagged protein-coding
#     pairs 0.795 true against 0.627 for the other family pairs (its own bar 0.90 not passed); non-coding pairs are unscored.
#   tools/rustle_pipeline.sh trunc5 --bam B --fasta G --out PREFIX [--threads N]   OPT-IN (register 1677): PREFIX.trunc5.gtf = PREFIX.gtf without
#                              multi-exon chains that at least as many primary reads run through and continue 5' of as end with
#   tools/rustle_pipeline.sh chainscore --bam B --fasta G --out PREFIX [--threads N]   (register 1761; applied to PREFIX.final.gtf by default since 2026-11-16, RUSTLE_CHAINSCORE=0 = off): PREFIX.chainscore.gtf = PREFIX.gtf without
#                              multi-exon chains whose frozen chain score is below the fixed cut (tools/data/chainscore_r1769.json, cut 0.02 since 2026-11-18, register 1810; RUSTLE_CHAINSCORE_TAU=0.055 = the earlier cut; default model since 2026-11-18)
#                              RUSTLE_CHAINSCORE_MODEL=r1677 = the earlier weight file (register 1677, cut 0.061), also for assemble's filter on PREFIX.final.gtf
#   tools/rustle_pipeline.sh capstarts --bam B --fasta G --out PREFIX [--threads N]   (register 1717 / 1720; also a last step of `all`, RUSTLE_CAPSTARTS=0 = off; writes nothing on a library without cap G): PREFIX.capstarts.gtf = PREFIX.gtf plus
#                              nested starts backed by >= 3 capped reads (non-templated 5' G); refuses libraries without a usable cap G
#   tools/rustle_pipeline.sh shiftdrop --bam B --fasta G --out PREFIX [--threads N]   OPT-IN (register 1728-1730 / 1732 / 1734 / 1738 / 1739): PREFIX.shiftdrop.gtf = PREFIX.gtf without
#                              weaker copies of a stronger chain with one internal exon end moved whose reads carry an indel beside the moved
#                              intron (>= half of them); 0.3-0.9% of chains, 2-4% of them annotated; no other stage reads this file
#   RUSTLE_SHORTEXONS (row 1878, user 2026-11-22): `assemble` named on its own repairs reads that skip an exon of <= 15 bp before it
#                              assembles (unset = 1; inside `all` unset = 0; 1 / 0 force it); needs python3 with pysam and room for a
#                              second BAM copy while it runs. PREFIX.ranked.gtf (rows 1880 / 1884; written by the oneread step,
#                              RUSTLE_RANKED_LIST=0 = off) = one ranked list: PREFIX.final.gtf by reads, then the one-read chains by score;
#                              RUSTLE_RANKED_TAIL=1 (opt-in) appends the chains the dropping rules remove (one more assembly)
#   tools/rustle_pipeline.sh shortexons --bam B --fasta G --out PREFIX [--threads N]   OPT-IN (register 1753 / 1754 / 1757): PREFIX.shortexons.gtf = the default
#                              assembly of the reads after those that skip an exon of <= 15 bp are rewritten through it (another read spells
#                              the exon, or >= 2 reads carry its bases and one splice-flanked site of the intron fits); repairs a small class,
#                              no consistent effect on precision; needs python3 with pysam; PREFIX.gtf is not touched. The inner assemble is the
#                              default one (RUSTLE_* switches apply, command-line assemble options do not); while it runs, PREFIX.shortexons.d
#                              holds a second copy of the BAM (peak about twice its size), removed at the end
#   tools/rustle_pipeline.sh earlyends --bam B --fasta G --out PREFIX [--threads N]   (register 1721 / 1722; also the last step of `all`, RUSTLE_EARLYENDS=0 = off): PREFIX.earlyends.gtf = PREFIX.gtf plus
#                              nested early ends whose reads end at one position with a polyA signal and no A-rich run
#   tools/rustle_pipeline.sh cache-ls --out PREFIX      list what PREFIX.cache holds (cache-clear: delete it)
#   tools/rustle_pipeline.sh sd --fasta G --out PREFIX [--threads N] [--bin DIR] [--no-cache]   the genome's SDs, cached (stage sd)
# --piecewise (catalog only; needs the cache): the catalog's representatives (both BAM passes + the span-overlap
#   collapse and its POA, none of which compares two contigs) are built one CONTIG per piece, each cached in
#   PREFIX.cache/reps/; when all are cached they are merged into the one-run order and the catalog continues as one
#   run (the k11 all-vs-all of ALL representatives, edges, families; families still cross contigs). Same products
#   as without it (cmp-checked 2026-09-25). --max-pieces N / --budget-s S bound one call: it computes that much,
#   appends to PREFIX.catalog.log and EXITS 75 while work remains (pieces pending, or the merge still to do) — call
#   it again until it exits 0. --piece-records R also cuts a contig of more than R BAM records into pieces of ~R
#   at positions no record crosses (still exact; human chr1 alone exceeds 10 min as one piece); --piece LABEL
#   computes that one piece only (as the log labels it, e.g. chr13:65273-15760749), for a known-heavy piece. The all-vs-all
#   itself can be sharded through RUSTLE_MINIMAP2 (tools/mm2_shard.sh).
# CACHE (default on, --no-cache turns it off): families and catalog keep their expensive intermediates in
#   PREFIX.cache/ (RUSTLE_CACHE_DIR): the catalog's collapsed representatives (reps/<key>/reps.tsv + reps.fa, both
#   BAM passes and the locus collapse) and every all-vs-all PAF (paf/<key>/out.paf). A re-run that changes only
#   downstream settings (the E_r edge rule, gamma, coverage split) replays them: human chr16 358 s -> 0.9 s,
#   byte-identical outputs. Keys cover the binary, BAM/FASTA (+ indexes) and every upstream RUSTLE_* setting.
#   The families PAF is keyed on every byte of PREFIX.fam.loci.fa (hashed as it is written) and replayed as a HARD
#   LINK to PREFIX.fam.loci.paf, not a copy (chimp: the hit step 2.6-13 s -> 0.01 s); a write through that link
#   invalidates the entry (mtime/inode/sampled-content pins), and RUSTLE_CACHE_VERIFY=1 re-hashes it on every hit.
# --inspect: also write the analyst dumps — catalog edge tables (PREFIX.cache/inspect/catalog.*: reps.fa, PAF,
#   edges.tsv, nodes.tsv, rule.tsv, params.tsv), per-round collapse statistics in the catalog log, and the
#   assignment evidence (PREFIX.assign.psv_*.tsv, PREFIX.assign.posterior.tsv).
# Loci are seeded from primaries PLUS secondaries within 2% of the molecule's GENOME-WIDE best alignment score
#   (one `as_table` pass over the BAM -> BAM.as_table.tsv next to the BAM, shared by every run on it, + its binary load
#   sidecar .asbin; PREFIX.molecules.tsv links to it; PREFIX.molecules.tsv itself when the BAM's directory is read-only): on gorilla NC_073244.2 (the
#   seeding pre-registration's verdict contig) this finds 30 more loci (97% annotated) and joins more
#   >=90%-identity referee pairs (21 -> 26 of 83 on genes with >= 2 exonic primary reads; all of the gain is one
#   tandem array) at pair precision 1.000 (rows 1060/1100/1101/1116); the streaming assembler applies the filter
#   at no memory cost (validated identical to the buffered path). --no-seed-secondaries restores primaries-only seeding (2026-09-24 default flip).
#   The table is only as genome-wide as the BAM: on a region SLICE it admits secondaries whose real best lies
#   outside the slice (tes44: 4,210 transcripts vs 3,911 with the full-BAM table) — give the driver the full BAM.
# SEEDING POOL, the switch (2026-10-07): which alignments seed the assembly. One flag, three pools, everything else equal:
#   --seed-pool primary   primary alignments only                                   (= --no-seed-secondaries)
#   --seed-pool good      primaries + secondaries with AS >= R x the molecule's genome-wide best AS, R = 0.98
#                         (the default of every stage but `assemble` on its own, which seeds with primary since register 1812;
#                         = --seed-secondaries; needs the as_table pass above)
#   RUSTLE_SEED_SPLICED=support (opt-in, passed through to the binary; docs/GOAL_multicopy_families_2026-11-19.md question 3): a
#                         SPLICED chain seeds only when >= 1 primary record carries it; secondaries still add support (streaming
#                         --assemble-only path only; the binary refuses it elsewhere). Default off = byte-identical.
#   RUSTLE_SEED_UNSPLICED=primary|support (opt-in, passed through to the binary; r1518/r1520): what the pool does with an UNSPLICED
#                         secondary: primary = drop it; support = it adds support to the primary-seeded cluster it overlaps but never
#                         seeds or bridges one (array-wide single-exon blobs gone; 13% of CAT single-exon pseudogene hits lost).
#   --seed-pool all       primaries + EVERY secondary (RUSTLE_GTF_SECONDARY=1, no filter, no table): the 2026-09-21 pool; in the
#                         2026-10-07 study 1.4-1.5x the transcripts of good on every contig, and the aligner
#                         records of good x17.5 (chr16; x22.9 of primary), x14.0 (chr17; x15.0), x3.1 (gorilla NPIP; x3.9)
#                         (human chr16: 14,183 against 9,473 transcripts, 1,648,588 against 94,433 PAF records;
#                         docs/archive/2026-11/LEDGER_O1.md:2696 (SEED_POOL_REAL_READS_2026-10-07.md))
#   --seed-as-ratio R     the width R of the good pool, a number in (0, 1]; 1 = exact ties only. Refused with another pool.
#   RUSTLE_SEED_POOL / RUSTLE_SEED_AS_RATIO set the same two as defaults (a flag wins; a malformed exported ratio is refused,
#   a valid one is ignored by the other pools). The pool the flag names is the pool the binary gets: an exported
#   RUSTLE_GTF_SECONDARY* of the shell is overridden (primary: unset; all: ratio and table unset), so those three variables
#   cannot contaminate a comparison between pools (other exported RUSTLE_* read-pool knobs still reach the binary;
#   bench/seed_pool/run.sh refuses any exported RUSTLE_*). One behaviour change: `--no-seed-secondaries` with an exported
#   RUSTLE_GTF_SECONDARY* used to let it through, and now unsets it (the default and the legacy flags are otherwise
#   byte-identical to the earlier driver, stub-checked). If the best-AS table is empty or from another BAM the binary says
#   so in PREFIX.assemble.log and the good pool fails with status 3 instead of running a wider pool. The
#   pool is announced in the log ("seeding pool: ..."), recorded in PREFIX.merged.env when it is not the default, and
#   changes nothing but that environment: the arguments of copy_assign are the same for every pool. Use a different
#   --out PREFIX per pool (bench/seed_pool/run.sh does, and prints the comparison). --as-table FILE reads (or builds, when absent
#   or made from another BAM) the best-AS table there instead of PREFIX.molecules.tsv, so many pools and widths share one pass over
#   --guide ANNOTATION.gtf[.gz] (assemble; opt-in, register 1622): also write PREFIX.guided.gtf = the de novo assembly with the
#   guide's transcripts that have a full-chain primary alignment added, and without de novo chains that are fragments or
#   intron-retaining versions of guide chains, have a new intron beside guide introns and < 10 reads, or lie outside the guide
#   with < 3 reads (`utilities guide`, src/guide.rs). PREFIX.gtf and every later stage stay de novo. Since register 1934
#   (2026-11-24) the rule reads PREFIX.final.gtf (PREFIX.gtf when there is none); RUSTLE_GUIDED_INPUT=raw = PREFIX.gtf as before.
#   the BAM. --contig NAME[,NAME...] (assemble) assembles those contigs only, copy_assign --region NAME:0-LEN (--regions FILE for
#   several) with LEN from FASTA.fai, in place of --genome-wide: the products equal the genome-wide run's on that contig (docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md) and it
#   makes a pool comparison on one contig a matter of minutes; the families, assign and flag stages read whatever assembly is there.
# The splice index is needed by `candidates` (the consensus sequences' genome hits) and `flag` (home search); the
# annotation (--gff) by `flag` only (IG/TR screen).
# Environment: RUSTLE_POLISH_SUBCHAIN=tag|drop adds `--polish-subchain` to `assemble` (default unset = off);
#   RUSTLE_POLISH_TSS=off|tag|rescue|split sets `--polish-tss` for `assemble` (default unset = rescue since 2026-11-01, register 1485);
#   RUSTLE_POLISH_TSS_END=median adds `--polish-tss-end` to `assemble` (default unset = off; register 1487);
#   RUSTLE_POLISH_TES_EXTEND=on adds `--polish-tes-extend` to `assemble` (default unset = off; register 1488);
#   RUSTLE_MONO_MODEL=cluster adds `--mono-model` to `assemble` (default unset = off; register 1489);
#   RUSTLE_FINAL_SUBCHAIN=drop|off: the sub-chain drop as a fourth final.gtf layer (default drop, register 1498);
#   RUSTLE_PYTHON=python3: the interpreter for families' extra representatives (RUSTLE_FAMREPS, default on) and assign's skip lists
#     (RUSTLE_EXCLUDE_BED / RUSTLE_EXCLUDE_SOFTMASK / RUSTLE_EXCLUDE_FASTA / RUSTLE_ANNOTATED_ONLY); not needed otherwise.
#   RUSTLE_EXCLUDE_SOFTMASK=0.9 (opt-in, register 1690): assign skips families whose copies are mostly soft-masked (>= that
#     share of exon bases lowercase in --fasta): satellite and transposon-class families, on any masked genome, no BED needed.
#     Does not catch rRNA reliably. May be set together with RUSTLE_EXCLUDE_BED (the lists are joined).
#   RUSTLE_EXCLUDE_BED=BED (rows 1559 / 1749 / 1750): assign skips families whose copies mostly lie in the BED's intervals (a
#     satellite-only list; docs/DATA.md). Unset: FASTA.exclude.bed when that sidecar exists (CHM13 has one since 2026-11-15, user);
#     empty: off.
#   RUSTLE_SOLE_CANDIDATE=1 (opt-in): assign reports a contested molecule with one candidate copy as `assigned`
#     (sole_candidate = 1), the default before 2026-11-13. Default now: such a molecule is `tied` (user ruling, registers
#     787 / 1656 / 1695); on A119b that is about 40% of the earlier assigned rows.
#   RUSTLE_EXCLUDE_FASTA=REF.fa (opt-in, register 1691): assign skips families whose copies are mostly covered (>= 50% of the
#     copy sequence) by minimap2 asm20 alignments to REF.fa. For rRNA / rDNA families, which the soft mask does not reach:
#     tools/data/rdna_human.fa (one human 45S unit with its spacer + one 5S unit; finds human and gorilla rDNA families) or
#     tools/data/rdna_human.masked.fa (the unit's own repeats masked: lists no other family, misses Alu-rich spacer transcripts).
#     The lists of the three RUSTLE_EXCLUDE_* switches are joined.
#   RUSTLE_ASSEMBLY_LAYERS=on|off: `assemble` also writes PREFIX.final.gtf, the layered assembly, from a second
#   parallel run (default on; PREFIX.gtf / PREFIX.families.gtf unchanged by construction; registers 1487-1491);
#   RUSTLE_FINAL_END_LAYERS="tssend tesext monomodel" (or `all`; default unset = none since row 1846, 2026-11-20; `all` = the
#     output before): the end layers of PREFIX.final.gtf only: median 5' ends (register 1487), 3' UTR extension (1488),
#     single-exon models from end clusters (1489). Row 1831: each helps only where it was tuned. The sub-chain layer
#     (RUSTLE_FINAL_SUBCHAIN) is not one of them and stays on;
#   RUSTLE_POLISH_TES=tag|pas-end adds `--polish-tes` to `assemble` (default unset = off);
#   RUSTLE_POLISH_MONO_PRIMED=0 adds `--no-polish-mono-primed` to `assemble` (default unset = 1 = on since 2026-10-26);
#   RUSTLE_POLISH_JUNCTION_SNAP=off|equiv|reads sets `--polish-junction-snap` of `assemble` (default unset = reads since
#     2026-11-18, row 1800; needs the BAM index);
#   RUSTLE_POLISH_SPLICED_SHARE=X sets `--polish-spliced-share X` of `assemble` (default unset = 0.05 since row 1848,
#     2026-11-21; 0 = off, the command before; r1801, row 1847): drops a
#     transcript with an intron spliced by less than X of the primary reads across it (0.05 in the r1801 diagnosis);
#   RUSTLE_POLISH_MONO_SHADOW_RATIO=R adds `--polish-mono-shadow-ratio R` (opt-in, r1801; unset or 0 = off): the shadow rule
#     keeps a single-exon transcript with at least R times the reads of the spliced transcript shadowing it;
#   RUSTLE_POLISH_ISM_RATIO_3P=R adds `--polish-ism-ratio-3p R` (opt-in, row 1843; unset = off): read ratio of the sub-chain
#     collapse for chains cut at the 3' side only (no longer chain reaches upstream of them); 0 = spared
#   RUSTLE_POLISH_JUNCTION_RATIO=R sets `--polish-junction-ratio R` (default unset = 20 since row 1845, 2026-11-20; 0 = off,
#     the command before): after the polish an intron is re-placed by the neighbour within 10 bp that has at least R times
#     its support (rows 1819, 1835);
#   RUSTLE_GTF_REGROUP=1 adds `--gtf-regroup` to `assemble` (RG3 regroup after polish; default unset = off; needs
#     RUSTLE_BRIDGE_REGROUP=off, whose default contains RG3's split);
#   RUSTLE_BRIDGE_REGROUP=off|f1|f1v2|f1units sets `--bridge-regroup` on `assemble` (unset = f1v2, THE DEFAULT since
#     2026-09-29; off = the 2026-09-25 products, byte for byte). With f1, f1v2 or f1units, `families`, and `flag` without
#     --gff, read PREFIX.families.gtf (f1 / f1v2: the assembly without its bridge transcripts, bridges are relations, not
#     loci; f1units, OPT-IN: each bridge transcript replaced by its UNITS, the exon runs between its bridge introns, see
#     copy_assign --help); every stage derives the mode the same way (BRIDGE_MODE below);
#   RUSTLE_BRIDGE_UNITS_LIST=FILE (f1units only; no whitespace in the path) adds `--bridge-units-list FILE` to `assemble`:
#     the transcripts and introns to cut come from FILE (any detector, e.g. bench/units_from_annotation.py) and F1's own
#     read evidence is not consulted;
#   RUSTLE_FAMILY_RELATIONS=1 adds `--emit-relations` to `families` (PREFIX.fam.relations.tsv, PREFIX.fam.members_by_locus.tsv:
#     one relation record per split transcript and the members of each family by locus; default unset = off);
#   RUSTLE_MIN_COV_SHORTER=C sets `mcl_families --min-cov-shorter C` on `families` (unset = the binary's own default,
#     0.70 since 2026-09-29, nothing passed; 0 = the 2026-09-25 edge weights);
#   RUSTLE_REPRESENTATIVE=most-reads|most-junctions sets `mcl_families --representative` on `families`, the transcript that
#     represents each de novo locus (loci.gff3, the copy table). Unset = nothing passed: the binary's default, most-reads
#     (the transcript with the most reads), byte for byte the earlier products; most-junctions = the transcript with the most
#     junctions (gaps >= 50 bp between exons), ties to the most reads — OPT-IN: the arm of
#     docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md, whose held-out decision (2026-10-04) kept most-reads: chr6 clause (a)
#     FAILED, Compara F 0.444 -> 0.372, on pre-f1v2 GTFs (docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md); an mcl_families
#     without the flag is refused;
#   RUSTLE_FAMILY_CONTAINER=1 adds `--emit-container` to `families` (PREFIX.fam.container*.tsv; default unset = off).
#   RUSTLE_FAMILY_EDGE_RULE=share|coverage (default share since 2026-11-09, row 1593): `families` builds edges by exon sharing
#     (`--min-shared-exon-frac 0.0 --min-exon-share ${RUSTLE_MIN_EXON_SHARE:-0.5} --inflation ${RUSTLE_INFLATION:-2.0}`);
#     coverage = the earlier rule (`--min-shared-exon-frac ${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}`, inflation 2.8), also chosen when
#     RUSTLE_MIN_SHARED_EXON_FRAC is set. The share rule aligns more pairs (exon gate at fraction 0): about 2x the families wall;
#   RUSTLE_SEDEF_BED=FILE turns on SD damping in `families` (unset: FASTA.sedef.bed when that sidecar exists; empty: off) (`--sd-damp-bed FILE --sd-damp ${RUSTLE_SD_DAMP:-0.25}`): edges
#     between loci no SEDEF pair links are down-weighted before MCL (row 1580: HG002 held-out +0.160 precision at equal
#     recall; A119b .802/.740 -> .925/.718). Set it whenever the genome has a SEDEF build (human CHM13, gorilla; paths in
#     docs/DATA.md); RUSTLE_SD_DAMP=1 keeps the BED but turns damping off. Unset with no FASTA.sedef.bed sidecar, or set empty = off, the earlier products.
# `families` is the DE NOVO mode (loci from the assembled GTF). The GUIDED mode (loci = the annotation's gene and
# pseudogene bodies, PREREG_heldout_families_2026-09-20 §2) is not a driver stage: figures/_o1_recovery.py
# (guided_families) runs its recipe step by step. Every product carries the PREFIX.
set -euo pipefail
# Run from a snapshot of this file (round-9 audit): bash reads a script as it executes it, so an edit or `git pull` of the
# repo during a multi-hour run spliced old and new code (seen 2026-11-09: "unexpected EOF" after a 10-minute assemble).
# The snapshot keeps size and mtime (cp -p), so the resume records' driver line is unchanged; RP_DRIVER_DIR (no RUSTLE_
# prefix: they must not enter the resume records) keeps the repo's tools/ directory for the sibling scripts and --bin.
# A stage call of `all` runs the same snapshot (its $0 is RP_SNAPSHOT); a call of the repo file (a user, rustle_samples.sh) snapshots.
if [ "${RP_SNAPSHOT:-}" != "$0" ]; then
  export RP_DRIVER_DIR RP_SNAPSHOT RP_OWNER_PID=$$; RP_DRIVER_DIR=$(cd "$(dirname "$0")" && pwd)
  RP_SNAPSHOT=$(mktemp "${TMPDIR:-/tmp}/rustle_pipeline.XXXXXX") && cp -p "$0" "$RP_SNAPSHOT" && exec bash "$RP_SNAPSHOT" "$@" \
    || { echo "[rustle_pipeline] cannot write the driver snapshot in ${TMPDIR:-/tmp}" >&2; exit 1; }
fi
# the process that made the snapshot (exec keeps its pid) removes it at exit; stage calls (other pids) leave it to their parent
if [ "${RP_OWNER_PID:-}" = "$$" ]; then trap 'rm -f "$RP_SNAPSHOT"' EXIT; fi
STAGE=${1:-all}; shift || true
ORIG_ARGS=("$@")
# Multi-sample runs (a comma in --bam): tools/rustle_samples.sh takes the whole call (pool, then split). A `case "$STAGE"` block
# at column 0, so figures/samples.py's driver_stage_code leaves it out of every stage's hash; one BAM never enters it.
case "$STAGE" in
  *) for ((i = 0; i + 1 < ${#ORIG_ARGS[@]}; i++)); do
       if [ "${ORIG_ARGS[$i]}" = --bam ] && [[ "${ORIG_ARGS[$((i + 1))]}" == *,* ]]; then
         exec bash "$RP_DRIVER_DIR/rustle_samples.sh" "$STAGE" "${ORIG_ARGS[@]}"
       fi
     done;;
esac
BAM=""; FASTA=""; OUT=""; INDEX=""; GFF=""; GUIDE=""; NAMES=""; THREADS=4; THREADS_SET=0; BIN="$RP_DRIVER_DIR/../target/release"; CONFIRM=(); FOREIGN=(); SEED_POOL=${RUSTLE_SEED_POOL:-}; SEED_RATIO=${RUSTLE_SEED_AS_RATIO:-0.98}; SEED_RATIO_FLAG=0; AS_TABLE=""; CONTIG=""; CONTIG_FLAG=0; CACHE=1; INSPECT=0
PIECEWISE=0; MAX_PIECES=0; BUDGET_S=0; PIECE_RECORDS=0; PIECE=""
LEGACY_CATALOG=0; CANDIDATES=0; DELTA=0.00958; CAND_MAX=1000
# `all --to STAGE` stops after STAGE (the escape hatch: `--to families` = the loci and their copy table before any read is
# assigned, 2026-10-30 NEXT TARGETS item 2); `--from STAGE` starts there (the earlier products must exist).
TO_STAGE=; FROM_STAGE=
INDEX_GB=8; LOW_MEM=0; NO_FLAG=0; RESUME=1; MAX_WINDOW_RECORDS=${RUSTLE_MAX_WINDOW_RECORDS:-0}
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAM=$2; shift 2;; --fasta) FASTA=$2; shift 2;; --out) OUT=$2; shift 2;; --index) INDEX=$2; shift 2;;
    --guide) GUIDE=$2; shift 2;;
    --names) NAMES=$2; shift 2;;
    --partners) # BAMs of OTHER libraries of the same species, comma-separated (= RUSTLE_ONEREAD_PARTNERS, colon-separated; row 1808)
      for _p in ${2//,/ }; do [ -s "$_p" ] || { echo "[rustle_pipeline] --partners: no such BAM '$_p'" >&2; exit 2; }; done
      export RUSTLE_ONEREAD_PARTNERS=${2//,/:}; shift 2;;
    --gff) GFF=$2; shift 2;; --threads) THREADS=$2; THREADS_SET=1; shift 2;; --bin) BIN=$2; shift 2;;
    --confirm) CONFIRM+=(--confirm "$2"); shift 2;; --foreign) FOREIGN+=(--foreign "$2"); shift 2;;
    --seed-pool) SEED_POOL=$2; shift 2;; --seed-as-ratio) SEED_RATIO=$2; SEED_RATIO_FLAG=1; shift 2;; --as-table) AS_TABLE=$2; shift 2;;
    --contig) CONTIG=$2; CONTIG_FLAG=1; shift 2;;
    --seed-secondaries) SEED_POOL=good; shift;; --no-seed-secondaries) SEED_POOL=primary; shift;;
    --cache) CACHE=1; shift;; --no-cache) CACHE=0; shift;; --inspect) INSPECT=1; shift;;
    --piecewise) PIECEWISE=1; shift;; --max-pieces) MAX_PIECES=$2; shift 2;; --budget-s) BUDGET_S=$2; shift 2;;
    --piece-records) PIECE_RECORDS=$2; shift 2;; --piece) PIECE=$2; shift 2;;
    --legacy-catalog) LEGACY_CATALOG=1; shift;; --candidates) CANDIDATES=1; shift;; --no-candidates) CANDIDATES=0; shift;;
    --to) TO_STAGE=$2; shift 2;; --from) FROM_STAGE=$2; shift 2;;
    --delta) DELTA=$2; shift 2;; --cand-max-reads) CAND_MAX=$2; shift 2;;
    --index-gb) INDEX_GB=$2; shift 2;; --low-mem) LOW_MEM=1; shift;; --no-flag) NO_FLAG=1; shift;; --no-resume) RESUME=0; shift;;
    --max-window-records) MAX_WINDOW_RECORDS=$2; shift 2;;
    *) echo "unknown argument $1" >&2; exit 2;;
  esac
done
# JOB_THREADS: without --threads, every CPU of the job for the steps whose products do not depend on the thread count
# (families' all-vs-all, as_table, the flag scan and verdict, the index build; round-37 audit). The settings records keep
# $THREADS, so a resumed run is not rerun for it.
JOB_THREADS=$THREADS; [ "$THREADS_SET" = 1 ] || JOB_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
# --candidates and --legacy-catalog name two different copy tables for assign: refused together, whatever the stage. Written
# as a `case "$STAGE"` block, as the final dispatch is, so figures/samples.py's driver_stage_code leaves it out of every
# stage's code hash: it refuses a call and changes no product.
case "$STAGE" in
  *) if [ "$CANDIDATES" = 1 ] && [ "$LEGACY_CATALOG" = 1 ]; then
       echo "[rustle_pipeline] --candidates and --legacy-catalog exclude each other (assign reads the families' copy table with its candidates, or the legacy catalog): drop one" >&2; exit 2
     fi;;
esac
# Register 1812 (user, 2026-11-19): the stage `assemble` on its own seeds with primary alignments unless a pool is named (flag,
# RUSTLE_SEED_POOL, or --seed-as-ratio / --as-table, which belong to good); every other stage, `all` included, keeps good.
if [ -z "$SEED_POOL" ]; then if [ "$STAGE" = assemble ] && [ "$SEED_RATIO_FLAG" = 0 ] && [ -z "$AS_TABLE" ]; then SEED_POOL=primary; else SEED_POOL=good; fi; fi
# The seeding pool (see the SEEDING POOL note above): a misspelling must not quietly run another pool. Same `case "$STAGE"`
# form as the check above.
case "$STAGE" in
  *) case "$SEED_POOL" in
       primary|good|all) ;;
       *) echo "[rustle_pipeline] --seed-pool (RUSTLE_SEED_POOL) must be primary, good or all (got '$SEED_POOL')" >&2; exit 2;;
     esac
     [[ "$SEED_RATIO" =~ ^(1(\.0+)?|0?\.[0-9]*[1-9][0-9]*)$ ]] || { echo "[rustle_pipeline] --seed-as-ratio (RUSTLE_SEED_AS_RATIO) must be a number in (0, 1] (got '$SEED_RATIO')" >&2; exit 2; }
     if [ "$SEED_RATIO_FLAG" = 1 ] && [ "$SEED_POOL" != good ]; then
       echo "[rustle_pipeline] --seed-as-ratio sets the width of the good pool and only the good pool has one (--seed-pool $SEED_POOL): drop one" >&2; exit 2
     fi
     if [ "$CONTIG_FLAG" = 1 ] && [ -z "$CONTIG" ]; then
       echo "[rustle_pipeline] --contig is empty (an empty shell variable?): name a contig, or drop the flag for a genome-wide run" >&2; exit 2
     fi;;
esac
if [ "$STAGE" = cache-clear ]; then
  [ -n "$OUT" ] || { echo "cache-clear needs --out PREFIX" >&2; exit 2; }
  [ -d "$OUT.cache" ] && du -sh "$OUT.cache" && rm -rf "$OUT.cache"
  exit 0
fi
if [ "$STAGE" = cache-ls ]; then
  [ -n "$OUT" ] || { echo "cache-ls needs --out PREFIX" >&2; exit 2; }
  for d in "$OUT".cache/*/*/; do
    [ -f "$d/DONE" ] || continue
    printf '%s\t%s\t%s\n' "$(du -sh "$d" | cut -f1)" "${d%/}" "$(head -1 "$d/key.tsv")"
  done
  exit 0
fi
# sd (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md): the genome's segmental duplications with the
# vendored sedef-rs (`utilities sd-scan`), once per genome: cached in RUSTLE_SD_CACHE_DIR (default PREFIX.cache) under
# sd/<md5 of the FASTA>-<md5 of sd-scan's version key, 8 hex>/genome.sd.*. PREFIX.sd.prefix names the cached prefix
# that `mcl_families --emit-hierarchy --sd-cache` takes; it is removed first and written only after a hit or a finished
# scan, so a failed run leaves none. Needs --fasta and --out only (no BAM), so it is a stage function plus a
# `case "$STAGE"` block that exits before the shared --bam check: figures/samples.py's driver_stage_code leaves every
# other stage's code hash as it was.
# The entry is marked done by DONE, written last: one `name<TAB>bytes` line per file of the entry (every genome.sd.* file
# and key.tsv), as src/family.rs's cache entries do. A hit needs the cache on (--no-cache always scans, into the same
# entry), a non-empty DONE and every listed file present with its listed size; anything else is a miss (an empty DONE,
# as this stage wrote before 2026-10-14, included). A miss removes DONE, scans, writes key.tsv, syncs the files, writes
# DONE.tmp and renames it to DONE.
stage_sd() {
  [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "[rustle_pipeline] sd needs --fasta and --out" >&2; exit 2; }
  [ -s "$FASTA" ] || { echo "[rustle_pipeline] sd: $FASTA is missing or empty" >&2; exit 2; }
  rm -f "$OUT.sd.prefix"
  local key md5 dir hit=0 f="" n="" b=""
  key=$("$BIN/utilities" sd-scan --version-key) || { echo "[rustle_pipeline] sd: $BIN/utilities has no sd-scan; rebuild it" >&2; exit 2; }
  md5=$(md5sum < "$FASTA" | cut -c1-32)
  dir=${RUSTLE_SD_CACHE_DIR:-$OUT.cache}/sd/$md5-$(printf '%s' "$key" | md5sum | cut -c1-8)
  # one stage per cache entry at a time (held until this stage exits): two misses on a shared cache used to rescan into
  # the same files, truncating the tables the run that finished first was reading
  # (opened for append, read-only when another user's lock file is not ours to write: a shared-cache hit then still runs;
  # NFS needs the writable fd for an exclusive lock, a read-only one locks shared)
  mkdir -p "${dir%/*}"
  # (the braces keep the 2>/dev/null to this open: on a bare exec it would silence the rest of the stage)
  local ro=0
  if { exec 8>> "$dir.lock"; } 2>/dev/null; then flock 8; else exec 8< "$dir.lock" && flock -s 8; ro=1; fi
  if [ "$CACHE" = 1 ] && [ -s "$dir/DONE" ]; then
    hit=1
    while IFS=$'\t' read -r n b || [ -n "$n" ]; do
      [ -n "$n" ] && [ -f "$dir/$n" ] && [ "$(wc -c < "$dir/$n")" = "$b" ] || { hit=0; break; }
    done < "$dir/DONE"
  fi
  if [ "$hit" = 1 ]; then
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: cached ($dir)" >&2
  else
    # a shared lock (the lock file is not ours to write) does not keep a second scanner out: scan only under our own lock
    [ "$ro" = 0 ] || { echo "[rustle_pipeline] sd: $dir is not a cache hit and $dir.lock is not writable by this user: scan with your own RUSTLE_SD_CACHE_DIR" >&2; exit 2; }
    mkdir -p "$dir"
    rm -f "$dir/DONE" "$dir/DONE.tmp"
    echo "[rustle_pipeline] $(date +%H:%M:%S) sd: utilities sd-scan on $FASTA -> $dir" >&2
    "$BIN/utilities" sd-scan --genome "$FASTA" --out "$dir/genome" --threads "$THREADS" > "$OUT.sd.log" 2>&1 \
      || { echo "[rustle_pipeline] sd: sd-scan failed (exit $?); see $OUT.sd.log" >&2; exit 1; }
    printf 'sd\t%s\t%s\t%s\n' "$md5" "$key" "$FASTA" > "$dir/key.tsv"
    sync "$dir"/genome.sd.* "$dir/key.tsv" || { echo "[rustle_pipeline] sd: could not sync the files of $dir; not marked done" >&2; exit 1; }
    for f in "$dir"/genome.sd.* "$dir/key.tsv"; do
      if [ -f "$f" ]; then printf '%s\t%s\n' "${f##*/}" "$(wc -c < "$f")"; fi
    done > "$dir/DONE.tmp"
    sync "$dir/DONE.tmp" && mv "$dir/DONE.tmp" "$dir/DONE" || { echo "[rustle_pipeline] sd: could not write $dir/DONE" >&2; exit 1; }
  fi
  printf '%s\n' "$dir/genome" > "$OUT.sd.prefix"
  echo "[rustle_pipeline] $(date +%H:%M:%S) sd: $(awk -F'\t' '!/^#/ && $11=="ok"' "$dir/genome.sd.pairs.tsv" | wc -l) SD pairs (ok), $(awk '!/^#/' "$dir/genome.sd.families.tsv" | wc -l) duplication families; --sd-cache $dir/genome ($OUT.sd.prefix)" >&2
}
case "$STAGE" in
  sd) stage_sd; exit 0;;
esac
if [ "$STAGE" = index ]; then
  [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "index needs --fasta and --out" >&2; exit 2; }
elif [ "$STAGE" = names ]; then
  [ -n "$OUT" ] || { echo "names needs --names ANNOT and --out" >&2; exit 2; }
else
  [ -n "$BAM" ] && [ -n "$FASTA" ] && [ -n "$OUT" ] || { echo "need --bam, --fasta, --out" >&2; exit 2; }
fi
# One writer per PREFIX (audit 2026-11-18): two sessions running `all --out X` deleted each other's X.final.* and both wrote
# X.gtf. Only the invocation the user started takes the lock; the stages it runs as children inherit it. A `case "$STAGE"`
# block at column 0, so figures/samples.py's driver_stage_code leaves it out of every stage's hash.
# shortcut: a background index build of a failed run still holds the lock, so a rerun is refused until it ends.
case "$STAGE" in
  all|merged|assemble|families|assign) if [ "${RP_OWNER_PID:-}" = "$$" ]; then
      mkdir -p "$(dirname "$OUT")"
      exec 7>> "$OUT.lock"; flock -n 7 || { echo "[rustle_pipeline] another run is writing $OUT.* (lock $OUT.lock)" >&2; exit 2; }
    fi;;
esac
# --low-mem (docs/archive/2026-11/LEDGER_O1.md:3502 (GENOMEWIDE_FAST_2026-10-12.md)): the splice index in 1-Gbp parts (flag holds one part at a time) and one
# assign region at a time; MM_THREADS caps families' all-vs-all and flag's align at 16 threads (their minimap2 memory scales with
# threads, outputs identical; gorilla families 38 -> 18 GB, flag align 44 -> 25 GB; audit rounds 43, 45); every other setting unchanged
[[ "$MAX_WINDOW_RECORDS" =~ ^[0-9]+$ ]] || { echo "[rustle_pipeline] --max-window-records must be a whole number (got '$MAX_WINDOW_RECORDS')" >&2; exit 2; }
if [ "$LOW_MEM" = 1 ]; then INDEX_GB=1; export RUSTLE_REGION_THREADS=1; [ "$MAX_WINDOW_RECORDS" -gt 0 ] || MAX_WINDOW_RECORDS=500000; fi
MM_THREADS=$JOB_THREADS; [ "$LOW_MEM" = 1 ] && [ "$MM_THREADS" -gt 16 ] && MM_THREADS=16
[[ "$INDEX_GB" =~ ^[1-9][0-9]*$ ]] || { echo "[rustle_pipeline] --index-gb must be a positive integer (got '$INDEX_GB')" >&2; exit 2; }
# skip lists that `merged` and the --candidates assign do not apply: refused here, before assemble and families run (round-23 audit);
# the stage functions refuse them again. Same `case "$STAGE"` form as the checks above.
case "$STAGE" in
  *) if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}${RUSTLE_EXCLUDE_SOFTMASK:-}${RUSTLE_EXCLUDE_FASTA:-}" ] \
        || { [ "$STAGE" = merged ] && [ -z "${RUSTLE_EXCLUDE_BED+x}" ] && [ -s "$FASTA.exclude.bed" ]; }; then
       if [ "$STAGE" = merged ]; then
         echo "[rustle_pipeline] merged applies no --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED / RUSTLE_EXCLUDE_SOFTMASK / RUSTLE_EXCLUDE_FASTA: use \`all\` for them" >&2; exit 2
       elif [ "$CANDIDATES" = 1 ] && { [ "$STAGE" = all ] || [ "$STAGE" = assign ]; }; then
         echo "[rustle_pipeline] --candidates does not apply --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED / RUSTLE_EXCLUDE_SOFTMASK / RUSTLE_EXCLUDE_FASTA: drop one of them" >&2; exit 2
       fi
     fi;;
esac
export TMPDIR=${TMPDIR:-/tmp}
if [ "$CACHE" = 1 ]; then export RUSTLE_CACHE_DIR="$OUT.cache"; else unset RUSTLE_CACHE_DIR; fi
INSPECT_CAT=(); INSPECT_ASSIGN=()
# ASSIGN_FAST (docs/archive/2026-11/LEDGER_O1.md:3502 (GENOMEWIDE_FAST_2026-10-12.md)): every copy_assign of the assign stage skips the POA homology DIAGNOSTIC
# (byte-identical families / assignments / abundance; ~85% of the wall time on dense families; only the diagnostic
# .fallback.tsv and its log counts are lost) and runs RUSTLE_REGION_THREADS regions at once (byte-identical). Defaults since
# r1543 (docs/PREREG_pipeline_speed_2026-11-05.md gw7, A119b genome-wide: 6,589 s vs 19,358 s, peak 199 GB on a 250 GB node):
# 8 regions, 4 family threads x 16 star threads, 4M-record read budget (24 regions / 16M on a large node with the satellite
# exclusion: assign_parallel_defaults below). LOW_MEM=1 keeps regions sequential.
# RUSTLE_POA_DIAGNOSTIC=1 restores the diagnostic pass.
# Round 54 audit: assign's regions at once and read budget, "REGIONS BUDGET". 24 / 16M on the configuration measured on A119b
# (>= 48 cores, >= 128 GB, RUSTLE_EXCLUDE_BED set: 286 s against 484 s, tables identical, peak RSS 40 GB against 34 GB; 36 / 32M
# gains nothing more). Elsewhere r1543's 8 / 4M: without the exclusion those peaked at 199 GB on the satellite families.
assign_parallel_defaults() {
  # CORES MEM_KB EXCLUDE_BED LOW_MEM
  if [ "$1" -ge 48 ] && [ "$2" -ge 134217728 ] && [ -n "$3" ] && [ "$4" != 1 ]; then echo "24 16000000"; else echo "8 4000000"; fi
}
# Round 58 audit: the same for assemble's sweep, "CONTIGS BUDGET": 24 / 80M on a large node (A119b 291 s at 8 / 20M, 237 s at
# 16 / 40M, 212 s at 25 / 80M; the three GTFs identical; 9 / 11 / 13 GB per process, two processes), else 8 / 20M.
assemble_parallel_defaults() {
  # CORES MEM_KB LOW_MEM
  if [ "$1" -ge 48 ] && [ "$2" -ge 134217728 ] && [ "$3" != 1 ]; then echo "24 80000000"; else echo "8 20000000"; fi
}
# Round 60: the binaries size a reader's inflate workers from the machine's cores (bam::bgzf_workers); under a smaller
# --threads / Slurm allocation the job's threads are the measure. Same rule: a quarter, at least min(4, threads), at most 16.
if [ -z "${RUSTLE_BAM_THREADS:-}" ]; then
  bw=$(( JOB_THREADS / 4 )); [ "$bw" -le 16 ] || bw=16; [ "$bw" -ge 4 ] || bw=$(( JOB_THREADS < 4 ? JOB_THREADS : 4 ))
  export RUSTLE_BAM_THREADS=$bw
fi
read -r ASM_REGIONS ASM_BUDGET < <(assemble_parallel_defaults "$JOB_THREADS" "$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)" "$LOW_MEM")
case "${RUSTLE_REGION_THREADS:-8}" in ''|*[!0-9]*|0) echo "[rustle_pipeline] RUSTLE_REGION_THREADS must be a positive integer (got '${RUSTLE_REGION_THREADS:-}')" >&2; exit 2;; esac
# RUSTLE_SOLE_CANDIDATE=1: assign with copy_assign --sole-candidate (the default before 2026-11-13)
SOLE=(); [ "${RUSTLE_SOLE_CANDIDATE:-0}" = 1 ] && SOLE=(--sole-candidate)
read -r ASSIGN_REGIONS ASSIGN_BUDGET < <(assign_parallel_defaults "$JOB_THREADS" "$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)" "${RUSTLE_EXCLUDE_BED:-}${RUSTLE_EXCLUDE_SOFTMASK:-}" "$LOW_MEM")
ASSIGN_FAST=(--region-threads "${RUSTLE_REGION_THREADS:-$ASSIGN_REGIONS}")
export RUSTLE_FAMILY_THREADS="${RUSTLE_FAMILY_THREADS:-4}" RUSTLE_STAR_THREADS="${RUSTLE_STAR_THREADS:-16}" RUSTLE_REGION_READ_BUDGET="${RUSTLE_REGION_READ_BUDGET:-$ASSIGN_BUDGET}"
[ "${RUSTLE_POA_DIAGNOSTIC:-0}" = 1 ] || ASSIGN_FAST+=(--skip-poa-diagnostic)
if [ "$INSPECT" = 1 ]; then
  mkdir -p "$OUT.cache/inspect"
  INSPECT_CAT=("RUSTLE_ER_EDGE_DUMP=$OUT.cache/inspect/catalog" RUSTLE_COLLAPSE_STATS=1)
  INSPECT_ASSIGN=(--dump-psv --posterior)
fi
# RUSTLE_POLISH_ISOFORM_FRACTION (unset = 0.02, the value this command always had): the isoform fraction of the polish. A switch only so
# that the relaxed assembly of RUSTLE_RANKED_TAIL=1 can turn the rule off (0); no other value has been judged as a default (rows 1837, 1876).
case "${RUSTLE_POLISH_ISOFORM_FRACTION:-0.02}" in ""|*[!0-9.]*|*.*.*) echo "[rustle_pipeline] RUSTLE_POLISH_ISOFORM_FRACTION must be a number >= 0 (got '${RUSTLE_POLISH_ISOFORM_FRACTION:-}')" >&2; exit 2 ;; esac
# RUSTLE_ASSEMBLY_POLISH=none (default full): the transcripts before the polishing rules, the input of tools/polish_rules.py
# (the Python version of those rules, for reading and for the reasons in the provenance tables); never for a real run.
# The assembler refuses junction snapping, the 5' rescue and junction replacement without its polishing (they run around
# it), so `none` switches those three off: the file is the transcripts BEFORE snapping, not the exact input of the rules in a
# normal run (there snapping comes first). Exact per-transcript reasons of a normal run would need the assembler to write them.
if [ "${RUSTLE_ASSEMBLY_POLISH:-full}" = none ]; then RUSTLE_POLISH_JUNCTION_SNAP=off; RUSTLE_POLISH_JUNCTION_RATIO=0; RUSTLE_POLISH_TSS=off; fi
case "${RUSTLE_ASSEMBLY_POLISH:-full}" in full|none) ;; *) echo "[rustle_pipeline] RUSTLE_ASSEMBLY_POLISH must be full or none (got '$RUSTLE_ASSEMBLY_POLISH')" >&2; exit 2 ;; esac
POLISH="--assembly-polish ${RUSTLE_ASSEMBLY_POLISH:-full} --polish-isoform-fraction ${RUSTLE_POLISH_ISOFORM_FRACTION:-0.02} --polish-mono-shadow --polish-mono-quantile ${RUSTLE_POLISH_MONO_QUANTILE:-0.82} --polish-ism-ratio ${RUSTLE_POLISH_ISM_RATIO:-0.3} --polish-retained-ratio ${RUSTLE_POLISH_RETAINED_RATIO:-10}"
# RUSTLE_POLISH_ISM_RATIO=R (register 1810; unset = 0.3 since 2026-11-18, 0.7 = the command before): a transcript whose intron chain is a
# contiguous part of a longer assembled chain is kept when it has at least R times that chain's reads. 0.7 was overfit to the primate
# samples (register 1793: it cost 9-10% of the correct chains on mouse); 0.1 is the non-primate choice of register 1800.
# RUSTLE_POLISH_RETAINED_RATIO=R (register 1745; unset = 10, the same command): a transcript whose exon contains an intron of its
# gene is dropped when that intron has >= R x the transcript's reads; 0 turns the retained-intron filter off.
# RUSTLE_POLISH_SUBCHAIN=tag|drop (opt-in; unset = off, the same command): copy_assign --polish-subchain, which tags
# (or drops) each transcript that is an end-compatible contiguous sub-chain of a longer emitted transcript of the same
# locus with >= 1/2 its reads (`subchain_of` / `subchain_missing`; drop is a documented dev-only trade, see its --help)
case "${RUSTLE_POLISH_SUBCHAIN:-}" in
  "") ;;
  off|tag|drop) POLISH="$POLISH --polish-subchain $RUSTLE_POLISH_SUBCHAIN" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_SUBCHAIN must be off, tag or drop (got '$RUSTLE_POLISH_SUBCHAIN')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_NASCENT=tag|drop (opt-in, register 1635; unset = off, the same command): copy_assign --polish-nascent, StringTie's
# nascent link (a transcript that follows a better-supported one's first introns, then runs into its next intron)
case "${RUSTLE_POLISH_NASCENT:-}" in
  "") ;;
  off|tag|drop) POLISH="$POLISH --polish-nascent $RUSTLE_POLISH_NASCENT" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_NASCENT must be off, tag or drop (got '$RUSTLE_POLISH_NASCENT')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TSS=off|tag|rescue|split: copy_assign --polish-tss, the read-proven TSS (tag: `tss_clusters` only; rescue:
# also keeps proven short-TSS forms the polish dropped; split: also emits a chain once per proven TSS cluster). Unset = rescue,
# the default since 2026-11-01 (register 1485, PREREG ASSEMBLY ENDS E1: held-out A119b chr11 / chr7 chains +7 / +9 at
# transcript precision -0.07 / -0.01 points); RUSTLE_POLISH_TSS=off restores the earlier output (`--polish-tss off` is byte-identical to no flag).
POLISH_TSS=${RUSTLE_POLISH_TSS:-rescue}
case "$POLISH_TSS" in
  off|tag|rescue|split) POLISH="$POLISH --polish-tss $POLISH_TSS" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TSS must be off, tag, rescue or split (got '$RUSTLE_POLISH_TSS')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TSS_END=median (opt-in; unset = off, the same command): copy_assign --polish-tss-end, each multi-exon transcript's
# 5' end at the lower median of its exact-chain reads' starts (register 1487, PASSED held-out chr11 / chr7: starts within 50 bp
# of CAT 40.6 -> 64.6% / 39.3 -> 63.4%, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_POLISH_TSS_END:-}" in
  "") ;;
  off|median) POLISH="$POLISH --polish-tss-end $RUSTLE_POLISH_TSS_END" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TSS_END must be off or median (got '$RUSTLE_POLISH_TSS_END')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TES_EXTEND=on (opt-in; unset = off, the same command): copy_assign --polish-tes-extend, a multi-exon model's 3' end
# extended to a PAS-proven cluster of unspliced reads over >= 2-read coverage (register 1488, PASSED held-out chr11 / chr7: base
# sensitivity +2.1 / +2.4, precision unchanged, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_POLISH_TES_EXTEND:-}" in
  "") ;;
  off|on) POLISH="$POLISH --polish-tes-extend $RUSTLE_POLISH_TES_EXTEND" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TES_EXTEND must be off or on (got '$RUSTLE_POLISH_TES_EXTEND')" >&2; exit 2 ;;
esac
# RUSTLE_MONO_MODEL=cluster (opt-in; unset = off, the same command): copy_assign --mono-model, single-exon models rebuilt from
# end-clustered unspliced reads (register 1489, PASSED held-out chr11 / chr7: base precision +17.2 / +16.5, transcript precision
# +0.2 / +0.4, chains identical); needs --polish-tss (the default rescue). Default pending the user.
case "${RUSTLE_MONO_MODEL:-}" in
  "") ;;
  off|cluster) POLISH="$POLISH --mono-model $RUSTLE_MONO_MODEL" ;;
  *) echo "[rustle_pipeline] RUSTLE_MONO_MODEL must be off or cluster (got '$RUSTLE_MONO_MODEL')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_TES=tag|pas-end (opt-in; unset = off, the same command): copy_assign --polish-tes, the read + genome
# proven TES (tag: `tes_clusters` / `tes_pas` / `tes_primed` only; pas-end: also moves an internally primed 3' end to
# the most-3' PAS-proven cluster of the transcript's own reads). Dev-only, in-sample evidence: see its --help
case "${RUSTLE_POLISH_TES:-}" in
  "") ;;
  off|tag|pas-end) POLISH="$POLISH --polish-tes $RUSTLE_POLISH_TES" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_TES must be off, tag or pas-end (got '$RUSTLE_POLISH_TES')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_MONO_PRIMED (default ON since 2026-10-26, user decision; unset or 1 = on, the same command; 0 = off):
# copy_assign drops a single-exon transcript when >= 80% of its unspliced reads are internally primed (StringTie's
# long-read poly-A filter, genomic form; held-out PASS bench/stgaps/PREREG_mono_primed_2026-10-26.md)
case "${RUSTLE_POLISH_MONO_PRIMED:-}" in
  ""|1) ;;
  0) POLISH="$POLISH --no-polish-mono-primed" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_MONO_PRIMED must be 0 or 1 (got '$RUSTLE_POLISH_MONO_PRIMED')" >&2; exit 2 ;;
esac
# RUSTLE_COMPAT_SUPPORT=1 (opt-in; unset or 0 = off, the same command): copy_assign --compat-support, compatible-read support in
# the assembler's chain admission (a chain short of the read floor is admitted when every junction is carried by enough reads
# whose chains are contiguous, end-compatible sub-chains of it; docs/archive/2026-11/LEDGER_O1.md:3055 (PREREG_assembler_compatible_support_2026-10-11.md))
case "${RUSTLE_COMPAT_SUPPORT:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --compat-support" ;;
  *) echo "[rustle_pipeline] RUSTLE_COMPAT_SUPPORT must be 0 or 1 (got '$RUSTLE_COMPAT_SUPPORT')" >&2; exit 2 ;;
esac
# RUSTLE_DIRECT_TRANSFRAGS=1|union (opt-in; unset or 0 = off, the same command): copy_assign --direct-transfrags (or, with union,
# --direct-transfrags-union: plus every chain the default admits), transcripts from
# maximal transfrags with the reads of the transfrags they contain transferred to them, no flow
# (docs/archive/2026-11/LEDGER_O1.md:3131 (PREREG_direct_transfrags_2026-10-11.md)); exclusive with RUSTLE_COMPAT_SUPPORT
case "${RUSTLE_DIRECT_TRANSFRAGS:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --direct-transfrags" ;;
  union) POLISH="$POLISH --direct-transfrags-union" ;;
  *) echo "[rustle_pipeline] RUSTLE_DIRECT_TRANSFRAGS must be 0, 1 or union (got '$RUSTLE_DIRECT_TRANSFRAGS')" >&2; exit 2 ;;
esac
# Opt-in levers toward FLAIR's exclusive transcripts (docs/archive/2026-11/LEDGER_GENERAL.md:511 (PENDING_2026-10-11.md); unset = off, the same command):
#   RUSTLE_TRANSFER_CREDIT=full|best copy_assign --transfer-credit full|best (with RUSTLE_DIRECT_TRANSFRAGS)
#   RUSTLE_READ_JUNCTION_SNAP=W    copy_assign --read-junction-snap W (bp; read junctions snap to better-supported ones)
#   RUSTLE_MIN_CHAIN_READS=N       copy_assign --min-chain-reads N (1 = single-read spliced chains)
#   RUSTLE_READ_STRAND=0           assemble: single-exon strand back to the '+' placeholder (default 1 = read-orientation vote, r1525)
#   RUSTLE_DUP_ADMIT=1             pass 1: coordinate-duplicate spliced records count toward the floor for a chain that would otherwise
#                                  fall under it (r1524: rescues 2-read family chains; genome-wide too broad, 68,749 chains at 1% CAT match)
#   RUSTLE_DUP_COUNT=0|1           pass 1 of `assemble` (default unset = 1 since row 1850, 2026-11-21; 0 = before; rows 1836, 1849):
#                                  coordinate-duplicate spliced records count in `reads` (isoform fraction, sub-chain ratio) but
#                                  not toward the 2-read floor, which stays on distinct records. The floor-1 run of the
#                                  `oneread` stage keeps 0 unless the variable is set
#   RUSTLE_CHECK_SPLICE=1          copy_assign --check-splice (FLAIR --check_splice on each spliced read)
#   RUSTLE_TERMINAL_EXON_ANCHOR=N  copy_assign --terminal-exon-anchor N (bp; FLAIR --stringent's end rule is 25)
#   RUSTLE_EXTEND_TRANSFRAGS=1     OPT-IN again since 2026-11-19 (user; register 1811: primate-only, +2.5% / +3.6% correct chains on
#                                  OR6737 / A119b, under 0.4% on five non-primate datasets with 7 or more wrong chains per correct
#                                  one). Default on from 2026-11-16 to 2026-11-19 (register 1760: its new chains beat the default's 2-read chains
#                                  on four unseen cell-line samples, +580 to +1,352 correct chains for -0.001 to -0.004 precision):
#                                  copy_assign --extend-transfrags (StringTie-style greedy extension, no flow); unset or 0 = off
#   RUSTLE_TRANSFER_ENDS=1         copy_assign --transfer-ends (ends from reads sharing the first / last junctions)
#   RUSTLE_READ_ISOFORM_K=K        copy_assign --read-isoform-k K (also admit a chain whose every junction >= K reads carry;
#                                  2 = StringTie's single-read choice, register 1267; buffered reader, slower)
case "${RUSTLE_TRANSFER_CREDIT:-}" in ""|split) ;; full|best) POLISH="$POLISH --transfer-credit $RUSTLE_TRANSFER_CREDIT" ;;
  *) echo "[rustle_pipeline] RUSTLE_TRANSFER_CREDIT must be split, full or best (got '$RUSTLE_TRANSFER_CREDIT')" >&2; exit 2 ;; esac
case "${RUSTLE_READ_JUNCTION_SNAP:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_READ_JUNCTION_SNAP must be a whole number of bp (got '$RUSTLE_READ_JUNCTION_SNAP')" >&2; exit 2 ;;
  *) POLISH="$POLISH --read-junction-snap $RUSTLE_READ_JUNCTION_SNAP" ;; esac
case "${RUSTLE_MIN_CHAIN_READS:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_MIN_CHAIN_READS must be a whole number (got '$RUSTLE_MIN_CHAIN_READS')" >&2; exit 2 ;;
  *) POLISH="$POLISH --min-chain-reads $RUSTLE_MIN_CHAIN_READS" ;; esac
case "${RUSTLE_CHECK_SPLICE:-}" in ""|0) ;; 1) POLISH="$POLISH --check-splice" ;;
  *) echo "[rustle_pipeline] RUSTLE_CHECK_SPLICE must be 0 or 1 (got '$RUSTLE_CHECK_SPLICE')" >&2; exit 2 ;; esac
case "${RUSTLE_TERMINAL_EXON_ANCHOR:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_TERMINAL_EXON_ANCHOR must be a whole number of bp (got '$RUSTLE_TERMINAL_EXON_ANCHOR')" >&2; exit 2 ;;
  *) POLISH="$POLISH --terminal-exon-anchor $RUSTLE_TERMINAL_EXON_ANCHOR" ;; esac
case "${RUSTLE_EXTEND_TRANSFRAGS:-}" in ""|0) ;; 1) POLISH="$POLISH --extend-transfrags" ;;
  *) echo "[rustle_pipeline] RUSTLE_EXTEND_TRANSFRAGS must be 0 or 1 (got '$RUSTLE_EXTEND_TRANSFRAGS')" >&2; exit 2 ;; esac
case "${RUSTLE_TRANSFER_ENDS:-}" in ""|0) ;; 1) POLISH="$POLISH --transfer-ends" ;;
  *) echo "[rustle_pipeline] RUSTLE_TRANSFER_ENDS must be 0 or 1 (got '$RUSTLE_TRANSFER_ENDS')" >&2; exit 2 ;; esac
case "${RUSTLE_READ_ISOFORM_K:-}" in ""|0) ;; *[!0-9]*) echo "[rustle_pipeline] RUSTLE_READ_ISOFORM_K must be a whole number (got '$RUSTLE_READ_ISOFORM_K')" >&2; exit 2 ;;
  *) POLISH="$POLISH --read-isoform-k $RUSTLE_READ_ISOFORM_K" ;; esac
# RUSTLE_POLISH_JUNCTION_SNAP=off|equiv|reads (unset = reads, the default since row 1800; off = the command before it): copy_assign --polish-junction-snap,
# the fuzzy-junction merge that keeps read-proven splice sites (reads: a junction within 10 bp of a better-supported one
# in its locus moves onto it only when its own reads do not prove it and the other's reads do). Evidence: rows 1793 and
# 1800 (docs/PREREG_snap_ism_defaults_2026-11-18.md); the chain score filter was fitted without it
case "${RUSTLE_POLISH_JUNCTION_SNAP:=reads}" in
  off|equiv|reads) POLISH="$POLISH --polish-junction-snap $RUSTLE_POLISH_JUNCTION_SNAP" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_JUNCTION_SNAP must be off, equiv or reads (got '$RUSTLE_POLISH_JUNCTION_SNAP')" >&2; exit 2 ;;
esac
# RUSTLE_POLISH_SPLICED_SHARE / RUSTLE_POLISH_MONO_SHADOW_RATIO (opt-in, r1801; unset or 0 = off, the same command)
# RUSTLE_POLISH_SPLICED_SHARE: unset = 0.05, the default since row 1848 (user 2026-11-21, after row 1847); 0 = off, the command before
case "${RUSTLE_POLISH_SPLICED_SHARE:=0.05}" in 0) ;; *[!0-9.]*) echo "[rustle_pipeline] RUSTLE_POLISH_SPLICED_SHARE must be a number >= 0 (got '$RUSTLE_POLISH_SPLICED_SHARE')" >&2; exit 2 ;;
  *) POLISH="$POLISH --polish-spliced-share $RUSTLE_POLISH_SPLICED_SHARE" ;; esac
case "${RUSTLE_POLISH_MONO_SHADOW_RATIO:-}" in ""|0) ;; *[!0-9.]*) echo "[rustle_pipeline] RUSTLE_POLISH_MONO_SHADOW_RATIO must be a number >= 0 (got '$RUSTLE_POLISH_MONO_SHADOW_RATIO')" >&2; exit 2 ;;
  *) POLISH="$POLISH --polish-mono-shadow-ratio $RUSTLE_POLISH_MONO_SHADOW_RATIO" ;; esac
# RUSTLE_POLISH_MONO_SPARE=1 / RUSTLE_POLISH_MONO_PAIR=1 (opt-in, row 1897 / 1901, user 2026-11-22; unset or 0 = off, the same command):
# pairing of single-exon transcripts with the multi-exon ones around them. SPARE: the shadow rule drops a single-exon transcript
# only when it overlaps an exon of a multi-exon transcript of its own strand (it can be a piece of it); one wholly in an intron,
# or covered only by an opposite-strand gene, is kept (6-43% annotated on four real libraries against 1.5-3.7% for the pieces,
# below the single-exon transcripts the default keeps: a recall option). PAIR: multi-exon transcripts get frag3_models /
# frag3_reads = the dropped single-exon transcripts lying in their last exon and ending within 50 bp of their 3' end.
case "${RUSTLE_POLISH_MONO_SPARE:-0}" in 0) ;; 1) POLISH="$POLISH --polish-mono-shadow-spare" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_MONO_SPARE must be 0 or 1 (got '$RUSTLE_POLISH_MONO_SPARE')" >&2; exit 2 ;; esac
case "${RUSTLE_POLISH_MONO_PAIR:-0}" in 0) ;; 1) POLISH="$POLISH --polish-mono-pair" ;;
  *) echo "[rustle_pipeline] RUSTLE_POLISH_MONO_PAIR must be 0 or 1 (got '$RUSTLE_POLISH_MONO_PAIR')" >&2; exit 2 ;; esac
# RUSTLE_POLISH_ISM_RATIO_3P=R (opt-in, row 1843; unset = off, the same command): copy_assign --polish-ism-ratio-3p R, the read
# ratio of the sub-chain collapse for chains cut at the 3' side only (0 = spared)
case "${RUSTLE_POLISH_ISM_RATIO_3P:-}" in "") ;; *[!0-9.]*) echo "[rustle_pipeline] RUSTLE_POLISH_ISM_RATIO_3P must be a number >= 0 (got '$RUSTLE_POLISH_ISM_RATIO_3P')" >&2; exit 2 ;;
  *) POLISH="$POLISH --polish-ism-ratio-3p $RUSTLE_POLISH_ISM_RATIO_3P" ;; esac
# RUSTLE_POLISH_JUNCTION_RATIO=R (unset = 20, the default since row 1845, user 2026-11-20; 0 = off, the command before):
# copy_assign --polish-junction-ratio (rows 1819, 1835); the chain score filter was fitted without it
case "${RUSTLE_POLISH_JUNCTION_RATIO:=20}" in 0) ;; *[!0-9.]*) echo "[rustle_pipeline] RUSTLE_POLISH_JUNCTION_RATIO must be a number >= 1 (got '$RUSTLE_POLISH_JUNCTION_RATIO')" >&2; exit 2 ;;
  *) POLISH="$POLISH --polish-junction-ratio $RUSTLE_POLISH_JUNCTION_RATIO" ;; esac
# RUSTLE_GTF_REGROUP=1 (opt-in; unset or 0 = off, the same command): copy_assign --gtf-regroup, RG3 — after every polish
# step, split a gene_id whose surviving transcripts share no same-strand exonic base (a dropped readthrough bridge, or
# two pre-polish components that collided on one base tid); split-only, intron chains unchanged, the deeper piece keeps
# the name and the others become <gene_id>.rg<k>. Exclusive with the bridge regroup below, whose default (f1v2) contains
# this split: RUSTLE_GTF_REGROUP=1 needs RUSTLE_BRIDGE_REGROUP=off. See its --help
case "${RUSTLE_GTF_REGROUP:-}" in
  ""|0) ;;
  1) POLISH="$POLISH --gtf-regroup" ;;
  *) echo "[rustle_pipeline] RUSTLE_GTF_REGROUP must be 0 or 1 (got '$RUSTLE_GTF_REGROUP')" >&2; exit 2 ;;
esac
# RUSTLE_BRIDGE_REGROUP=off|f1|f1v2|f1units (unset = f1v2, THE DEFAULT since 2026-09-29 by the user's decision on the F1v2 held-out
# Outcome and its family-level side result, docs/archive/2026-09/PREREG_o1_cover_growth_2026-09-29.md; off = the 2026-09-25 products,
# byte for byte): copy_assign --bridge-regroup, passed explicitly by `assemble`. A transcript that is the only link
# between two pieces of its gene_id, when reads end at a PAS inside its intron and other reads start there at their own
# promoter (f1v2: and it carries fewer reads than each piece), becomes a relation `<gene_id>.fus<k>` with `fusion_of`,
# and the pieces split as RUSTLE_GTF_REGROUP splits them (so the two are exclusive). With f1 or f1v2, `assemble` also
# writes PREFIX.families.gtf, the GTF without the bridges, and `families` reads it, as the held-out runs did
# (docs/archive/2026-09/PREREG_f1_bridge_locus_2026-09-28.md, docs/archive/2026-09/PREREG_f1v2_readshare_2026-09-29.md), as does `flag`'s scan of the
# de novo loci (without --gff): a bridge is a relation, never a locus to scan. See its --help.
# BRIDGE_MODE is the one derivation every stage uses (assemble's flag, families' and flag's input, the guard below).
# f1units (OPT-IN, docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md Part C): F1's bridges without its read-share rule are
# kept in PREFIX.families.gtf as units instead of being left out (PREFIX.gtf is exactly f1's: only the transcripts that were cut are
# tagged); RUSTLE_BRIDGE_UNITS_LIST names the cuts instead of F1.
BRIDGE_MODE=${RUSTLE_BRIDGE_REGROUP:-f1v2}
FAM_GTF=$OUT.gtf
case "$BRIDGE_MODE" in
  off) ;;
  f1|f1v2|f1units)
    [ "${RUSTLE_GTF_REGROUP:-0}" = 0 ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP=$BRIDGE_MODE (unset = f1v2) already splits every gene_id as RUSTLE_GTF_REGROUP does: set RUSTLE_BRIDGE_REGROUP=off with RUSTLE_GTF_REGROUP=1" >&2; exit 2; }
    FAM_GTF=$OUT.families.gtf ;;
  *) echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP must be off, f1, f1v2 or f1units (got '$RUSTLE_BRIDGE_REGROUP')" >&2; exit 2 ;;
esac
POLISH="$POLISH --bridge-regroup $BRIDGE_MODE"
case "${RUSTLE_BRIDGE_UNITS_LIST:-}" in
  "") ;;
  *[[:space:]]*) echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST must not contain whitespace (got '$RUSTLE_BRIDGE_UNITS_LIST')" >&2; exit 2 ;;
  *) [ "$BRIDGE_MODE" = f1units ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST names the units of RUSTLE_BRIDGE_REGROUP=f1units (got '$BRIDGE_MODE')" >&2; exit 2; }
     [ -s "$RUSTLE_BRIDGE_UNITS_LIST" ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_UNITS_LIST=$RUSTLE_BRIDGE_UNITS_LIST is missing or empty" >&2; exit 2; }
     POLISH="$POLISH --bridge-units-list $RUSTLE_BRIDGE_UNITS_LIST" ;;
esac
# RUSTLE_MIN_COV_SHORTER=C (unset = nothing passed: the binary's default, 0.70 since 2026-09-29 by the user's decision,
# register 1006/1014; 0 = the 2026-09-25 edge weights, byte-identical to every catalog built before the flip):
# mcl_families --min-cov-shorter C, the §6x4 containment escape — a pair whose coverage of the LONGER locus fails also
# passes at coverage >= C of the SHORTER locus's exonic length, with that coverage as its edge weight. Known regressions
# (register 1007/1009): NPIP in GUIDED mode, Soto F .833 -> .800; semi-guided SD-region nodes, precision .973 -> .833.
FAM_COV=()
case "${RUSTLE_MIN_COV_SHORTER:-}" in
  "") ;;
  *) [[ "$RUSTLE_MIN_COV_SHORTER" =~ ^(0|1|1\.0+|0?\.[0-9]+)$ ]] || { echo "[rustle_pipeline] RUSTLE_MIN_COV_SHORTER must be a number in [0, 1] (got '$RUSTLE_MIN_COV_SHORTER')" >&2; exit 2; }
     FAM_COV=(--min-cov-shorter "$RUSTLE_MIN_COV_SHORTER") ;;
esac
# RUSTLE_FAMILY_CONTAINER=1 (opt-in; unset or 0 = off, the same command): mcl_families --emit-container, the container
# of each family member's extra pieces (docs/archive/2026-09/PREREG_fusion_container_sim_2026-09-28.md §1): every clustered locus's
# all-transcript exon blocks, core (aligned to an exon base of another member of its family) or accessory, and the
# other families each accessory block aligns to -> PREFIX.fam.container.tsv / .container_relations.tsv /
# .container_summary.tsv. It never changes a family. See its --help
FAM_EXTRA=()
case "${RUSTLE_FAMILY_CONTAINER:-}" in
  ""|0) ;;
  1) FAM_EXTRA=(--emit-container) ;;
  *) echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER must be 0 or 1 (got '$RUSTLE_FAMILY_CONTAINER')" >&2; exit 2 ;;
esac
# RUSTLE_FAMILY_RELATIONS=1 (opt-in; unset or 0 = off, the same command): mcl_families --emit-relations, the relation
# record of every unit-split transcript (RUSTLE_BRIDGE_REGROUP=f1units) and the members of each family by locus ->
# PREFIX.fam.relations.tsv / PREFIX.fam.members_by_locus.tsv (docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md Part C).
# It never changes a family. Without units the relations table is empty and every locus is a `whole` member.
FAM_REL=()
case "${RUSTLE_FAMILY_RELATIONS:-}" in
  ""|0) ;;
  1) FAM_REL=(--emit-relations) ;;
  *) echo "[rustle_pipeline] RUSTLE_FAMILY_RELATIONS must be 0 or 1 (got '$RUSTLE_FAMILY_RELATIONS')" >&2; exit 2 ;;
esac
say() { echo "[rustle_pipeline] $(date +%H:%M:%S) $*" >&2; }
# The de novo loci a stage reads ($FAM_GTF: PREFIX.families.gtf when BRIDGE_MODE is f1, f1v2 or f1units, else PREFIX.gtf) must
# come from the assembly on disk: with a bridge mode (the default), PREFIX.families.gtf must exist and be no older than
# PREFIX.gtf; with off, a newer PREFIX.families.gtf means the GTF was assembled with a bridge mode and its bridges would
# be read as loci. $1 = the stage.
fam_gtf_guard() {
  if [ "$FAM_GTF" != "$OUT.gtf" ]; then
    [ -s "$FAM_GTF" ] && [ ! "$FAM_GTF" -ot "$OUT.gtf" ] || { echo "[rustle_pipeline] RUSTLE_BRIDGE_REGROUP=$BRIDGE_MODE (unset = f1v2) but $FAM_GTF is missing or older than $OUT.gtf: run assemble with it, or set RUSTLE_BRIDGE_REGROUP=off for $1 on a GTF assembled with off" >&2; exit 2; }
  elif [ -e "$OUT.families.gtf" ] && [ ! "$OUT.families.gtf" -ot "$OUT.gtf" ]; then
    echo "[rustle_pipeline] $OUT.gtf was assembled with a bridge mode (its bridges are relations, not loci) but RUSTLE_BRIDGE_REGROUP=off: unset it, or set assemble's mode, for $1 too" >&2; exit 2
  fi
}

stage_assemble() {
  # RUSTLE_DUP_COUNT (row 1850, user 2026-11-21 after row 1849): unset = 1 for this stage, every coordinate-duplicate spliced
  # record counts in a chain's `reads` (the 2-read floor stays on distinct records); 0 = the collapsed counts before
  : "${RP_DUP_COUNT_USER=${RUSTLE_DUP_COUNT:-}}"  # what the caller set, for the floor-1 run of the oneread stage (it is called inside this scope too)
  local -x RUSTLE_DUP_COUNT="${RUSTLE_DUP_COUNT:-1}"
  case "$RUSTLE_DUP_COUNT" in 0|1) ;; *) echo "[rustle_pipeline] RUSTLE_DUP_COUNT must be 0 or 1 (got '$RUSTLE_DUP_COUNT')" >&2; exit 2 ;; esac
  # RUSTLE_SHORTEXONS (row 1878, user 2026-11-22): the short-exon repair is part of this stage. minimap2 drops exons of <= 15 bp
  # (row 1746); the reads that skip one are rewritten through it (bench/famcmp/short_exons.py, no annotation) and the assembly
  # runs once, on the rewritten reads. Unset = 1 when `assemble` is named, 0 inside `all` (families / assign read the original
  # BAM and were not measured on repaired models); 1 / 0 force it. Row 1878: correct chains and precision both up on 9 of 11
  # datasets (3 simulated, 4 mouse, chimpanzee, human), the two gorillas within 4 chains. Needs python3 with pysam and room
  # for a second copy of the BAM in $OUT.shortexons.d while the stage runs (RUSTLE_SHORTEXONS_KEEP=1 keeps it); side files
  # $OUT.shortexons.units.tsv / .reads.tsv / .log. A failed repair stops the stage: RUSTLE_SHORTEXONS=0 assembles without it.
  local se="${RUSTLE_SHORTEXONS:-}"; if [ -z "$se" ]; then if [ -n "${RP_ALL:-}" ]; then se=0; else se=1; fi; fi
  case "$se" in 0|1) ;; *) echo "[rustle_pipeline] RUSTLE_SHORTEXONS must be 0 or 1 (got '$se')" >&2; exit 2 ;; esac
  if [ "$se" = 1 ] && [ -z "${RP_SHORTEXONS_FROM:-}" ]; then
    local sd="$OUT.shortexons.d" rc=0; mkdir -p "$sd"
    say "assemble: short-exon repair of $BAM (RUSTLE_SHORTEXONS=0 = off)"
    ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/short_exons.py" "$BAM" "$FASTA" "$sd/reads.bam" "$OUT.shortexons" "$THREADS" > "$OUT.shortexons.log" 2>&1 \
      || { [ -n "${RUSTLE_SHORTEXONS_KEEP:-}" ] || rm -rf "${sd:?}"; echo "[rustle_pipeline] assemble: the short-exon repair failed (see $OUT.shortexons.log); RUSTLE_SHORTEXONS=0 assembles without it" >&2; exit 1; }
    say "assemble: $(tail -1 "$OUT.shortexons.log")"
    ( export RP_SHORTEXONS_FROM="$BAM"; BAM="$sd/reads.bam"; stage_assemble ) || rc=$?
    [ -n "${RUSTLE_SHORTEXONS_KEEP:-}" ] || rm -rf "${sd:?}"
    [ "$rc" = 0 ] || exit "$rc"
    return 0
  fi
  say "assemble: $BAM -> $OUT.gtf"
  # Pre-flight (warnings only, never fatal): untrimmed polyA tails or unoriented reads in the first 100,000 primary reads.
  "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/check_reads.py" "$BAM" "$FASTA" >&2 || say "assemble: read pre-flight check did not run (needs python3 + pysam); continuing"
  # --contig NAME[,NAME...]: these contigs only (copy_assign --region NAME:0-LEN, or --regions FILE for several, LEN from the FASTA
  # index), in place of --genome-wide; the products equal the genome-wide run's on those contigs (docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md gates)
  local scope=(--genome-wide)
  if [ -n "$CONTIG" ]; then
    case "$CONTIG" in *$'\n'*|*$'\r'*) echo "[rustle_pipeline] --contig: a newline or carriage return in the value (the list is comma separated)" >&2; exit 2;; esac
    case ",$CONTIG," in *,,*) echo "[rustle_pipeline] --contig '$CONTIG': an empty name in the list" >&2; exit 2;; esac
    local items=() c clen regions="" nctg=0 seen=","
    IFS=, read -ra items <<< "$CONTIG"
    for c in "${items[@]}"; do
      case "$c" in
        *:*|*[[:space:]]*) echo "[rustle_pipeline] --contig '$c': a name with ':' or whitespace cannot be a region (copy_assign splits at the first ':')" >&2; exit 2;;
      esac
      case "$seen" in *",$c,"*) echo "[rustle_pipeline] --contig lists $c twice" >&2; exit 2;; esac
      seen+="$c,"
      clen=""
      if [ -s "$FASTA.fai" ]; then clen=$(awk -F'\t' -v c="$c" '($1 "")==(c ""){print $2; exit}' "$FASTA.fai"); fi   # as strings: awk would equate 01, 1.0 and 1e0
      [ -n "$clen" ] || { echo "[rustle_pipeline] --contig $c: no such contig in $FASTA.fai (samtools faidx writes it)" >&2; exit 2; }
      regions+="$c:0-$clen"$'\n'; nctg=$((nctg + 1))
    done
    if [ "$nctg" = 1 ]; then
      scope=(--region "${regions%$'\n'}")
      say "assemble: contig $CONTIG only (--region ${regions%$'\n'})"
    else
      printf '%s' "$regions" > "$OUT.assemble_regions.txt"
      scope=(--regions "$OUT.assemble_regions.txt")
      say "assemble: contigs $CONTIG only ($nctg regions in $OUT.assemble_regions.txt)"
    fi
  fi
  # the seeding pool (the SEEDING POOL note above): the three RUSTLE_GTF_SECONDARY* variables the binary reads are set, or
  # unset, for every pool, so a variable exported in the shell cannot change the pool the flag names
  local seed_env=(-u RUSTLE_GTF_SECONDARY -u RUSTLE_GTF_SECONDARY_AS_RATIO -u RUSTLE_GTF_SECONDARY_AS_TABLE)
  if [ "$SEED_POOL" = good ]; then
    say "assemble: seeding pool: good (AS >= $SEED_RATIO x the molecule's genome-wide best)"
    # The table depends on the BAM alone, so by default it is SHARED by every run on that BAM: BAM.as_table.tsv next to
    # the BAM (2026-11-07; the pass is ~200 s of I/O on a 96 GB BAM, register 1567), or PREFIX.molecules.tsv when the BAM's
    # directory is not writable. --as-table FILE names it instead. Reused when its header names this BAM and it is not
    # older than the BAM; otherwise rebuilt under a temporary name and renamed into place (a concurrent run never reads
    # a half-written table; rename keeps the inode and mtime the .asbin sidecar is keyed on).
    local table=${AS_TABLE:-}
    as_table_ok() { head -1 "$1" 2>/dev/null | grep -qF "bam=$(readlink -f "$BAM")$(printf '\t')bam_bytes=$(stat -Lc %s "$BAM")$(printf '\t')" && [ ! "$1" -ot "$BAM" ]; }
    if [ -z "$table" ]; then
      # next to the BAM when a valid table is there or the directory takes a new file (-w alone misses a full file quota);
      # otherwise PREFIX.molecules.tsv, which the check below then reuses like any other table
      local probe="$BAM.as_table.probe.$$"
      table="$OUT.molecules.tsv"
      if as_table_ok "$BAM.as_table.tsv" || { : > "$probe"; } 2>/dev/null; then table="$BAM.as_table.tsv"; fi
      rm -f "$probe"
    fi
    if as_table_ok "$table"; then
      say "assemble: genome-wide best-AS table: reusing $table"
    else   # absent, truncated, from another BAM, or older than the BAM
      say "assemble: genome-wide best-AS table -> $table (one pass over the BAM)"
      local tmp="$table.tmp.$$"
      if "$BIN/as_table" --bam "$BAM" --out "$tmp" --threads "$JOB_THREADS" > "$OUT.as_table.log" 2>&1 && [ -s "$tmp" ]; then
        [ -e "$tmp.asbin" ] && mv -f "$tmp.asbin" "$table.asbin"
        mv -f "$tmp" "$table"
      else
        rm -f "$tmp" "$tmp.tmp" "$tmp.asbin"
        echo "[rustle_pipeline] as_table failed (see $OUT.as_table.log)" >&2; exit 3
      fi
    fi
    # PREFIX.molecules.tsv stays the per-run name later tools read (tools/rustle_reassemble.sh): a link to the shared table
    if [ "$(readlink -f "$table")" != "$(readlink -f "$OUT.molecules.tsv")" ]; then
      ln -sfn "$(readlink -f "$table")" "$OUT.molecules.tsv"
      ln -sfn "$(readlink -f "$table").asbin" "$OUT.molecules.tsv.asbin"
    fi
    seed_env=(RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO="$SEED_RATIO" "RUSTLE_GTF_SECONDARY_AS_TABLE=$table")
  elif [ "$SEED_POOL" = all ]; then
    say "assemble: seeding pool: all (every secondary alignment, no filter)"
    seed_env=(-u RUSTLE_GTF_SECONDARY_AS_RATIO -u RUSTLE_GTF_SECONDARY_AS_TABLE RUSTLE_GTF_SECONDARY=1)
  else
    say "assemble: seeding pool: primary (primary alignments only)"
  fi
  # r1525 PASS (2026-11-04): single-exon models take their strand from the read-orientation vote (default on here;
  # RUSTLE_READ_STRAND=0 restores the '+' placeholder). A119b: spliced chains identical, CAT single-exon gene hits 1,465 -> 1,490,
  # complete truth families 204 -> 209, 3,031 placeholder-strand models now shadowed by their same-strand spliced model.
  # After the pool branches on purpose: they reassign seed_env wholesale (the first ship, b8d29f0d, was inert for the good pool).
  seed_env+=("RUSTLE_READ_STRAND=${RUSTLE_READ_STRAND:-1}")
  # the exact flags of this run, for tools/rustle_reassemble.sh's pass 2 (it must not re-derive them from the defaults)
  printf '%s\n' "$POLISH" > "$OUT.assemble.polish"
  # ... and its seeding environment and contig scope, one word per line (pass 2 must seed and scope as pass 1 did)
  printf '%s\n' "${seed_env[@]}" > "$OUT.assemble.seed_env"
  printf '%s\n' "${scope[@]}" > "$OUT.assemble.scope"
  # PREFIX.final.gtf (2026-11-02, registers 1487-1491): the layered assembly. The three end layers (median 5' ends, 3' UTR extension,
  # single-exon models from end clusters; opt-in since row 1846) each PASSED on held-out chromosomes, but each changed the families stage when it read
  # them (r1490 / r1491). So the layered assembly is a SECOND run, in parallel, written only as PREFIX.final.gtf: PREFIX.gtf and
  # PREFIX.families.gtf (every later stage's input) are byte-identical to a run without it. RUSTLE_ASSEMBLY_LAYERS=off skips it.
  local layers=()
  case "${RUSTLE_ASSEMBLY_LAYERS:-on}" in
    on) # a layer the environment already set (RUSTLE_POLISH_TSS_END / _TES_EXTEND / RUSTLE_MONO_MODEL) is in $POLISH for both
        # runs, and a second copy of the flag stops clap: the user's value wins (audit 2026-11-07, fourth pass)
        # RUSTLE_FINAL_END_LAYERS (row 1846, user 2026-11-20): the three end layers are opt-in, named here; unset = none
        # (before: all three). Row 1831: each helps only where it was tuned (mono-model and 3' extension on ape testis,
        # median 5' ends on A119b and mouse ES).
        local fl w; fl=${RUSTLE_FINAL_END_LAYERS:-}; fl=${fl//,/ }
        for w in $fl; do
          case "$w" in
            all) fl="tssend tesext monomodel" ;;
            tssend|tesext|monomodel|none) ;;
            *) echo "[rustle_pipeline] RUSTLE_FINAL_END_LAYERS takes tssend, tesext, monomodel, all or none (got '$w')" >&2; exit 2 ;;
          esac
        done
        for w in $fl; do
          case "$w" in
            tssend) case "$POLISH" in *--polish-tss-end*) ;; *) layers+=(--polish-tss-end median) ;; esac ;;
            tesext) case "$POLISH" in *--polish-tes-extend*) ;; *) layers+=(--polish-tes-extend on) ;; esac ;;
            monomodel) case "$POLISH" in *--mono-model*) ;; *) layers+=(--mono-model cluster) ;; esac ;;
          esac
        done
        # RUSTLE_FINAL_SUBCHAIN=drop|off (register 1498 PASS; unset = drop): --polish-subchain drop as a fourth layer of PREFIX.final.gtf only
        case "${RUSTLE_FINAL_SUBCHAIN:-drop}" in
          drop) case "$POLISH" in *--polish-subchain*) ;; *) layers+=(--polish-subchain drop) ;; esac ;;
          off) ;;
          *) echo "[rustle_pipeline] RUSTLE_FINAL_SUBCHAIN must be drop or off (got '$RUSTLE_FINAL_SUBCHAIN')" >&2; exit 2 ;;
        esac ;;
    off) ;;
    *) echo "[rustle_pipeline] RUSTLE_ASSEMBLY_LAYERS must be on or off (got '$RUSTLE_ASSEMBLY_LAYERS')" >&2; exit 2 ;;
  esac
  # the three end layers read the --polish-tss evidence pass: with RUSTLE_POLISH_TSS=off only the sub-chain layer is kept
  # (it needs no --polish-tss; round-23 audit), and with no layer left there is no layered assembly to build
  if [ "${#layers[@]}" -gt 0 ] && [ "$POLISH_TSS" = off ]; then
    local kept=() i
    for (( i=0; i<${#layers[@]}; i+=2 )); do [ "${layers[i]}" = --polish-subchain ] && kept+=("${layers[i]}" "${layers[i+1]}"); done
    if [ "${#kept[@]}" -gt 0 ]; then say "assemble: RUSTLE_POLISH_TSS=off, so $OUT.final.gtf has only ${kept[*]} (the end layers need --polish-tss)"
    else say "assemble: RUSTLE_POLISH_TSS=off, so no layered $OUT.final.gtf (its end layers need --polish-tss); RUSTLE_ASSEMBLY_LAYERS=off silences this"; fi
    layers=("${kept[@]}")
  fi
  # r1564 PASS (2026-11-07): contigs swept in parallel (--region-threads; low-mem mode sets 1 = serial), A119b assemble
  # 2,588 -> 1,685 s, six products byte-identical, RSS 16 GB. The read budget is 20M records here (80M on a large node, round 58), not assign's 4M:
  # at 4M chr1/chr2/chr13 (5-6M records each) run alone and the sweep stays near serial.
  local fin_pid=""
  # the side products of the previous assembly go too: a failed or switched-off step must not leave them beside the new GTF
  rm -f "$OUT.final.gtf" "$OUT.final.pid" "$OUT.final.rc" "$OUT".final_run.* "$OUT.final.named.gtf" "$OUT.named.gtf" "$OUT.chainscore.gtf" "$OUT.earlyends.gtf" "$OUT.capstarts.gtf"
  if [ "${#layers[@]}" -gt 0 ]; then
    # one subshell = the run, then its product under its final name (a rename, so PREFIX.final.gtf is whole or absent)
    ( env "${seed_env[@]}" RUSTLE_REGION_READ_BUDGET="$ASM_BUDGET" "$BIN/copy_assign" --assemble-only --region-threads "${RUSTLE_REGION_THREADS:-$ASM_REGIONS}" "${scope[@]}" --assembly-junctions strict $POLISH "${layers[@]}" --gtf-tpm \
      --bam "$BAM" --fasta "$FASTA" --out "$OUT.final_run" > "$OUT.final.log" 2>&1 \
        && [ -s "$OUT.final_run.gtf" ] && mv -f "$OUT.final_run.gtf" "$OUT.final.gtf"; rc=$?; rm -f "$OUT".final_run.*; echo "$rc" > "$OUT.final.rc" ) &
    fin_pid=$!
  fi
  # set -e ends the stage when the base run fails: the background layered run must not outlive it holding memory and threads
  local rc=0
  # PREFIX.polish_reasons.tsv (2026-11-24, user): every transcript the polishing removed, with the rule that removed it, written
  # by the assembler itself (exact, unlike tools/polish_rules.py on a file made with polishing off). Read by stage `provenance`.
  # RUSTLE_POLISH_REASONS=0 = no table; a binary built before the option is detected and simply gets no such argument.
  local reasons=(); rm -f "$OUT.polish_reasons.tsv"
  if [ "${RUSTLE_POLISH_REASONS:-1}" != 0 ] && [ "${RUSTLE_ASSEMBLY_POLISH:-full}" = full ] && [[ "$("$BIN/copy_assign" --help 2>/dev/null)" == *--polish-reasons* ]]; then   # not `| grep -q`: under pipefail an early-closed pipe (the help is 100 kB) reads as "no such option"
    reasons=(--polish-reasons "$OUT.polish_reasons.tsv"); fi
  env "${seed_env[@]}" RUSTLE_REGION_READ_BUDGET="$ASM_BUDGET" "$BIN/copy_assign" --assemble-only --region-threads "${RUSTLE_REGION_THREADS:-$ASM_REGIONS}" "${scope[@]}" --assembly-junctions strict $POLISH "${reasons[@]}" --gtf-tpm \
    --bam "$BAM" --fasta "$FASTA" --out "$OUT" > "$OUT.assemble.log" 2>&1 || rc=$?
  if [ "$rc" != 0 ]; then
    if [ -n "$fin_pid" ]; then pkill -P "$fin_pid" 2>/dev/null || true; kill "$fin_pid" 2>/dev/null || true; wait "$fin_pid" 2>/dev/null || true; rm -f "$OUT".final_run.*; fi
    echo "[rustle_pipeline] assemble failed (exit $rc, see $OUT.assemble.log)" >&2; exit "$rc"
  fi
  if [ -n "$fin_pid" ] && [ "${RP_LAYERS_BG:-0}" = 1 ]; then
    # round 66: under `all` the layered run (about 25 s longer on A119b) ends under families; `all` waits for this pid
    # pid, host and start time: a later run must not take another process that reuses the pid for this one (audit 2026-11-19)
    echo "$fin_pid $(hostname) $(ps -o lstart= -p "$fin_pid" 2>/dev/null | tr -s ' ' '_')" > "$OUT.final.pid"
    say "assemble: the layered assembly $OUT.final.gtf is finishing in the background"
  elif [ -n "$fin_pid" ]; then
    wait "$fin_pid" || true; layered_report
  fi
  # the good pool is only the good pool if the binary used the table: an empty table is ignored and the region-local rule applies, a table
  # from another BAM keeps every secondary it does not know; either way the run is wider than --seed-pool good says
  if [ "$SEED_POOL" = good ] && grep -qE 'table IGNORED|the table may come from a different BAM' "$OUT.assemble.log"; then
    echo "[rustle_pipeline] --seed-pool good: copy_assign did not apply the best-AS table ($(grep -m1 -E 'table IGNORED|the table may come from a different BAM' "$OUT.assemble.log" | cut -c1-300)); the assembly is not the good pool: delete $table and $table.asbin, then run again (an absent table is rebuilt; see $OUT.assemble.log)" >&2; exit 3
  fi
  # copy_assign writes PREFIX.families.gtf and PREFIX.gtf a few ms apart; on a filesystem with sub-second mtimes fam_gtf_guard's
  # `-ot` can then call it older (seen on NFS, 2026-10-29). Stamp it with PREFIX.gtf's mtime: same assembly, same age.
  [ "$FAM_GTF" = "$OUT.gtf" ] || [ ! -s "$FAM_GTF" ] || touch -r "$OUT.gtf" "$FAM_GTF"
  say "assemble: $(awk -F'\t' '$3=="transcript"' "$OUT.gtf" | wc -l) transcripts"
  if [ "$BRIDGE_MODE" = f1units ]; then
    say "assemble: $(grep -c 'fusion_of "' "$OUT.gtf") bridge transcripts cut into $(awk -F'\t' '$3=="transcript"' "$FAM_GTF" | grep -c 'fusion_unit "') units in the families input $FAM_GTF ($OUT.bridge_units.tsv)"
  elif [ "$FAM_GTF" != "$OUT.gtf" ]; then
    say "assemble: $(grep -c 'fusion_of "' "$OUT.gtf") bridge transcripts kept as fusion_of relations; families input $FAM_GTF"
  fi
}
# Register 1816 (user, 2026-11-19): the default model is r1813 (register 1813: nothing taken from the sample, so a chain's score does
# not move with the rest of the chain set; RUSTLE_CHAINSCORE_MODEL=r1769 = the model before), and the 5' truncation rule is on.
# RUSTLE_TRUNC5_RATIO=R (register 1815; unset = 3.22, 0 = off): the same pass also drops a chain when the reads
# that carry it and splice on past its 5' end outnumber the reads matching it exactly by more than R. Needs the chain score pass on.
# Register 1761, default since 2026-11-16 (user; RUSTLE_CHAINSCORE=0 = off): the layered assembly loses the multi-exon chains the
# frozen chain score puts below its fixed cut. The score is computed on $OUT.gtf (what the rule was tested on) and the same chain
# list is taken out of $OUT.final.gtf;
# $OUT.chainscore.gtf is $OUT.gtf filtered. $OUT.gtf and the families input are never changed. A failure (no python3 with numpy and
# pysam, say) is reported and leaves $OUT.final.gtf unfiltered.
# shortcut: one more BAM pass after every assemble (1-2 min on a 12M-read sample); fold the features into copy_assign if it matters
# RUSTLE_CHAINSCORE_MODEL picks the weight file of the chain score filter: unset or r1769 = the default since 2026-11-18 (user;
# register 1769: decile model fit on 11 samples of six species, no ORF feature, cut 0.02 since register 1810 (was 0.055 = at least 0.95 of the correct chains kept
# on each of the 12 samples; HG002 unseen keeps 0.978 at precision 0.280 against 0.970 at 0.259 for r1677); r1677 = the default
# until then (fit on OR6737 + A119b, cut 0.061); or a path to a model file.
export RUSTLE_TRUNC5_RATIO=${RUSTLE_TRUNC5_RATIO:-3.22}
chainscore_model() {
  local m="${RUSTLE_CHAINSCORE_MODEL:-r1813}"
  case "$m" in */*|*.json) printf '%s\n' "$m" ;; *) printf '%s\n' "$RP_DRIVER_DIR/data/chainscore_$m.json" ;; esac
}
final_chainscore() {
  [ "${RUSTLE_CHAINSCORE:-1}" != 0 ] || return 0
  if ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/chainscore.py" apply "$(chainscore_model)" "$OUT.gtf" "$BAM" "$FASTA" \
       "$OUT.chainscore.gtf.tmp" "$THREADS" "$OUT.final.gtf" "$OUT.final.gtf.tmp" > "$OUT.chainscore.log" 2>&1 && [ -s "$OUT.final.gtf.tmp" ]; then
    mv -f "$OUT.chainscore.gtf.tmp" "$OUT.chainscore.gtf"; mv -f "$OUT.final.gtf.tmp" "$OUT.final.gtf"
    say "assemble: chain score filter, $(grep -m1 'by the same chain list' "$OUT.chainscore.log" | sed 's/^[^:]*: //') taken out of $OUT.final.gtf ($OUT.chainscore.gtf: $(tail -1 "$OUT.chainscore.log"))"
  else
    rm -f "$OUT.chainscore.gtf.tmp" "$OUT.final.gtf.tmp"
    say "assemble: WARNING the chain score filter failed (see $OUT.chainscore.log); $OUT.final.gtf is the unfiltered layered assembly"
  fi
}
# oneread: OPT-IN (register 1766, rule S3b; RUSTLE_ONEREAD=1 also appends its models to $OUT.final.gtf at the end of assemble).
# $OUT.oneread.gtf = one-read chains the default lacks that pass tools/oneread.py: every junction also in another read, > 0.05 of the
# reads on each junction colinear, no longer read, ORF >= 900 nt with no early stop. Costs a second assemble at a 1-read floor
# ($OUT.oneread.d/) and one BAM pass. $OUT.gtf and the families input are untouched.
# $OUT.oneread.ranked.gtf (register 1806) = ALL one-read chains the default lacks, best oneread_score first (bench/famcmp/oneread_tier.py,
# model tools/data/oneread_r1899.json, RUSTLE_ONEREAD_MODEL=r1806 = the model before): a separate tier to browse, not a filter; one more BAM pass.
# One floor-1 assembly for both the one-read tiers and the multi-copy stage (register 1937, user 2026-11-24). The two stages
# asked for the same run but for the pooling-guard bypass; with RUSTLE_MULTICOPY_GATE now 1 by default they are the same
# command, and whichever stage runs second takes the first one's GTF. f1_settings = everything the run depends on that a
# stage chooses; a floor-1 GTF is reused only when its recorded settings equal the caller's and it is newer than $OUT.gtf.
# RUSTLE_F1_SHARE=0 = each stage assembles for itself, as before.
f1_settings() { printf 'gate=%s mono=%s dup=%s pool=%s ratio=%s contig=%s astable=%s bam=%s polish=%s\n' "$1" "$2" "${RP_DUP_COUNT_USER:-0}" "$SEED_POOL" "$SEED_RATIO" "${CONTIG:-}" "${AS_TABLE:-}" "$BAM" "$POLISH"; }
f1_shared() {  # f1_shared OTHER_DIR SETTINGS -> prints OTHER_DIR/f1.gtf when it can stand for this stage's run
  [ "${RUSTLE_F1_SHARE:-1}" != 0 ] && [ -s "$1/f1.gtf" ] && [ -s "$1/f1.settings" ] && [ "$1/f1.gtf" -nt "$OUT.gtf" ] && [ "$(cat "$1/f1.settings")" = "$2" ] && echo "$1/f1.gtf"
  return 0
}
stage_oneread() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] oneread needs $OUT.gtf (run assemble first), --bam and --fasta" >&2; return 2; }
  mkdir -p "$OUT.oneread.d"
  # RUSTLE_ONEREAD_GATE (unset = 1, the default since row 1868, user 2026-11-21; 0 = the command before): the floor-1 run also sets
  # RUSTLE_GATE_MIN_READS=1, so a single-read chain the gate's pooling guard marks as bridging two loci (row 1864: full-length reads of
  # long genes) is a candidate too. Row 1866: 65-99% more correct promoted chains at 0.01-0.03 lower precision; the stage runs 5-39% longer
  local gate1=(); case "${RUSTLE_ONEREAD_GATE:-1}" in 0) ;; 1) gate1=(RUSTLE_GATE_MIN_READS=1) ;;
    *) echo "[rustle_pipeline] RUSTLE_ONEREAD_GATE must be 0 or 1 (got '$RUSTLE_ONEREAD_GATE')" >&2; return 2 ;; esac
  local f1set f1have; f1set=$(f1_settings "${RUSTLE_ONEREAD_GATE:-1}" 0); f1have=$(f1_shared "$OUT.multicopy.d" "$f1set")
  { if [ -n "$f1have" ]; then
      ln -sf "$(cd "$(dirname "$f1have")" && pwd)/f1.gtf" "$OUT.oneread.d/f1.gtf"; printf '%s\n' "$f1set" > "$OUT.oneread.d/f1.settings"
      echo "[rustle_pipeline] oneread: floor-1 assembly taken from $f1have (the multi-copy stage made it with the same settings)"
    else
      rm -f "$OUT.oneread.d/f1.gtf" "$OUT.oneread.d/f1.settings"
      env "${gate1[@]}" RUSTLE_DUP_COUNT="${RP_DUP_COUNT_USER:-0}" RUSTLE_MIN_CHAIN_READS=1 RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_ONEREAD=0 RUSTLE_SEED_POOL="$SEED_POOL" RUSTLE_SEED_AS_RATIO="$SEED_RATIO" bash "$0" assemble --bam "$BAM" --fasta "$FASTA" \
        --out "$OUT.oneread.d/f1" --threads "$THREADS" ${BIN:+--bin "$BIN"} ${CONTIG:+--contig "$CONTIG"} ${AS_TABLE:+--as-table "$AS_TABLE"} && printf '%s\n' "$f1set" > "$OUT.oneread.d/f1.settings"
    fi && ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/oneread.py" "$OUT.gtf" "$OUT.oneread.d/f1.gtf" "$BAM" "$FASTA" "$OUT.oneread.gtf.tmp"; } > "$OUT.oneread.log" 2>&1 \
    || { echo "[rustle_pipeline] oneread failed (see $OUT.oneread.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.oneread.gtf.tmp"; return 1; }
  mv -f "$OUT.oneread.gtf.tmp" "$OUT.oneread.gtf"
  # register 1806: every one-read chain the default lacks, ranked (a ranking, not a filter; never appended to final.gtf)
  # RUSTLE_ONEREAD_MODEL (row 1900, user 2026-11-22): unset = r1899, the 23 features of row 1806 refitted on twenty datasets at the
  # defaults of 2026-11-22; r1806 = the ranking before (fitted 2026-11-18 on seven datasets under collapsed read counts). Row 1899:
  # with each species held out of its fit the refit was NOT adoptable by its bar (at least the old order on 12 of 17 datasets:
  # primates 11 of 13, mouse 1 of 4); adopted by the user against it. It changes the ORDER of the one-read tier only.
  local om="${RUSTLE_ONEREAD_MODEL:-r1899}"
  case "$om" in r1899|r1806) ;; *) echo "[rustle_pipeline] RUSTLE_ONEREAD_MODEL must be r1899 or r1806 (got '$om')" >&2; return 2 ;; esac
  if ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/oneread_tier.py" apply "$RP_DRIVER_DIR/data/oneread_$om.json" "$OUT.gtf" "$OUT.oneread.d/f1.gtf" \
       "$BAM" "$FASTA" "$OUT.oneread.ranked.gtf.tmp" "$THREADS" > "$OUT.oneread.ranked.log" 2>&1; then
    mv -f "$OUT.oneread.ranked.gtf.tmp" "$OUT.oneread.ranked.gtf"; say "oneread: $(tail -1 "$OUT.oneread.ranked.log") -> $OUT.oneread.ranked.gtf"
  else rm -f "$OUT.oneread.ranked.gtf.tmp"; say "oneread: WARNING the ranked tier failed (see $OUT.oneread.ranked.log); $OUT.oneread.gtf is unaffected"; fi
  # $OUT.ranked.gtf (rows 1880 / 1884, user 2026-11-22; RUSTLE_RANKED_LIST=0 = off): ONE ranked list of the assembly, tier 1 = the
  # transcripts of $OUT.final.gtf by reads, tier 2 = the ranked one-read chains it lacks by oneread_score (tools/ranked_list.py:
  # list_tier, list_rank, list_sensitive = tier 1 + the best 10% of tier 2). No transcript changes. Cut at another assembler's
  # output size the list held at least its annotated chains on 65 of 67 comparisons (20 datasets, row 1880); tier 2 is a recall tier.
  case "${RUSTLE_RANKED_LIST:-1}" in 0|1) ;; *) echo "[rustle_pipeline] RUSTLE_RANKED_LIST must be 0 or 1 (got '$RUSTLE_RANKED_LIST')" >&2; return 2 ;; esac
  rm -f "$OUT.ranked.gtf"
  if [ "${RUSTLE_RANKED_LIST:-1}" = 1 ] && [ -s "$OUT.oneread.ranked.gtf" ] && [ -s "$OUT.final.gtf" ]; then
    # RUSTLE_RANKED_TAIL=1 (OPT-IN, row 1884, user 2026-11-22; unset = 0): a third tier at the end of the list = the multi-exon chains
    # of one more assembly with a one-read floor and the dropping rules off (isoform fraction, sub-chain collapse, retained-intron,
    # spliced share) that neither earlier tier holds, by reads. It costs that assembly; on HG002 / A119b it is about a million
    # chains at 1.4-1.8% annotated, and it is what puts the list ahead of FLAIR on HG002. Never part of any other product.
    local tail=()
    case "${RUSTLE_RANKED_TAIL:-0}" in 0) ;;
      1) if env RUSTLE_MIN_CHAIN_READS=1 RUSTLE_GATE_MIN_READS=1 RUSTLE_POLISH_ISOFORM_FRACTION=0 RUSTLE_POLISH_ISM_RATIO=0 RUSTLE_POLISH_RETAINED_RATIO=0 RUSTLE_POLISH_SPLICED_SHARE=0 \
              RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_ONEREAD=0 RUSTLE_SEED_POOL="$SEED_POOL" RUSTLE_SEED_AS_RATIO="$SEED_RATIO" bash "$0" assemble --bam "$BAM" --fasta "$FASTA" \
              --out "$OUT.oneread.d/tail" --threads "$THREADS" ${BIN:+--bin "$BIN"} ${CONTIG:+--contig "$CONTIG"} ${AS_TABLE:+--as-table "$AS_TABLE"} > "$OUT.ranked.tail.log" 2>&1 \
            && [ -s "$OUT.oneread.d/tail.gtf" ]; then tail=(--tail "$OUT.oneread.d/tail.gtf")
         else say "oneread: WARNING the relaxed assembly for the third tier failed (see $OUT.ranked.tail.log); the list has two tiers"; fi ;;
      *) echo "[rustle_pipeline] RUSTLE_RANKED_TAIL must be 0 or 1 (got '$RUSTLE_RANKED_TAIL')" >&2; return 2 ;; esac
    if ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/ranked_list.py" "$OUT.final.gtf" "$OUT.oneread.ranked.gtf" "$OUT.ranked.gtf.tmp" "${tail[@]}" > "$OUT.ranked.log" 2>&1; then
      mv -f "$OUT.ranked.gtf.tmp" "$OUT.ranked.gtf"; say "oneread: $(tail -1 "$OUT.ranked.log") -> $OUT.ranked.gtf"
    else rm -f "$OUT.ranked.gtf.tmp"; say "oneread: WARNING the ranked list failed (see $OUT.ranked.log); the other products are unaffected"; fi
    rm -f "$OUT.oneread.d"/tail.*   # the relaxed assembly is read once
  fi
  # register 1808, OPT-IN: RUSTLE_ONEREAD_PARTNERS = colon-separated BAMs of OTHER libraries of the same species. Ranked chains in the
  # best 20% that a partner read also carries -> $OUT.oneread.promoted.gtf (oneread_promoted "1"). A RECALL tier: 0.22-0.30 precision on
  # WTC11, mouse adrenal and chimpanzee testis under the defaults of 2026-11-21 (row 1866), below the default's two-read chains there
  # (0.25-0.31); row 1808's "as precise as the two-read chains" held under the defaults of 2026-11-18 only. It needs a second library.
  rm -f "$OUT.oneread.promoted.gtf"
  if [ -n "${RUSTLE_ONEREAD_PARTNERS:-}" ] && [ -s "$OUT.oneread.ranked.gtf" ]; then
    local IFS=:; local partners=($RUSTLE_ONEREAD_PARTNERS); unset IFS
    if ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/oneread_tier.py" promote "$OUT.oneread.ranked.gtf" "$OUT.oneread.gtf" "$OUT.oneread.promoted.gtf.tmp" \
         "${partners[@]}" > "$OUT.oneread.promoted.log" 2>&1; then
      mv -f "$OUT.oneread.promoted.gtf.tmp" "$OUT.oneread.promoted.gtf"; say "oneread: $(tail -1 "$OUT.oneread.promoted.log") -> $OUT.oneread.promoted.gtf"
    else rm -f "$OUT.oneread.promoted.gtf.tmp"; say "oneread: WARNING the promotion failed (see $OUT.oneread.promoted.log)"; fi
  fi
  say "oneread: $(tail -1 "$OUT.oneread.log") -> $OUT.oneread.gtf"
}
# RUSTLE_ONEREAD=1: the strict tier (and the promoted chains, if any) appended to $OUT.final.gtf. Partner libraries alone
# (--partners / RUSTLE_ONEREAD_PARTNERS, user 2026-11-20): only the promoted chains are appended (a recall tier: 0.22-0.30 precision under
# today's defaults, row 1866; row 1808's "as precise as the default's two-read chains" held on 2026-11-18 only); RUSTLE_ONEREAD=0 keeps $OUT.final.gtf free of one-read models whatever else is set.
final_oneread() {
  local on=${RUSTLE_ONEREAD:-}
  [ "$on" = 1 ] || { [ "$on" != 0 ] && [ -n "${RUSTLE_ONEREAD_PARTNERS:-}" ]; } || return 0
  if stage_oneread; then
    if [ "$on" = 1 ]; then cat "$OUT.oneread.gtf" >> "$OUT.final.gtf"; say "assemble: one-read models appended to $OUT.final.gtf"; fi
    if [ -s "$OUT.oneread.promoted.gtf" ]; then cat "$OUT.oneread.promoted.gtf" >> "$OUT.final.gtf"; say "assemble: promoted one-read models (oneread_promoted \"1\") appended to $OUT.final.gtf"; fi
  else say "assemble: WARNING the one-read step failed; $OUT.final.gtf has no one-read models"; fi
}
# Register 1766 (user 2026-11-16): transcript lines of $OUT.final.gtf get orf_len and ptc attributes (tools/gtf_ptc.py; ptc "1" = the
# longest ORF stops > 50 nt before the last junction). A tag only, nothing is dropped; RUSTLE_PTC=0 = off; a failure leaves the file as it was.
final_ptc() {
  [ "${RUSTLE_PTC:-1}" != 0 ] || return 0
  if ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/gtf_ptc.py" "$OUT.final.gtf" "$FASTA" "$OUT.final.gtf.ptc.tmp" > "$OUT.ptc.log" 2>&1 && [ -s "$OUT.final.gtf.ptc.tmp" ]; then
    mv -f "$OUT.final.gtf.ptc.tmp" "$OUT.final.gtf"; say "assemble: $(tail -1 "$OUT.ptc.log") in $OUT.final.gtf"
  else
    rm -f "$OUT.final.gtf.ptc.tmp"; say "assemble: WARNING the ptc tag failed (see $OUT.ptc.log); $OUT.final.gtf has no ptc attribute"
  fi
}
layered_report() {
  if [ -s "$OUT.final.gtf" ]; then
    say "assemble: $(awk -F'\t' '$3=="transcript"' "$OUT.final.gtf" | wc -l) transcripts in the layered assembly $OUT.final.gtf (log $OUT.final.log)"
    final_chainscore
    final_oneread
    final_ptc
  else
    say "assemble: WARNING the layered assembly failed (see $OUT.final.log); $OUT.gtf and the families input are unaffected"
  fi
}
# guided: --guide GTF (register 1622, rule G1d; opt-in), a second product beside the de novo assembly. $OUT.gtf and every later
# stage stay de novo; $OUT.guided.gtf = de novo models the guide does not explain away + guide transcripts with a full-chain
# read. Its own stage (audit round 46): `all` calls it whenever the product is missing or older than its inputs, so a resumed
# run that skips assemble still writes it.
# multicopy: OPT-IN, named only (registers 1858-1867; never part of `all`). The multi-copy mode's second half; the first half is
# `assemble --seed-pool good`. Copies of a family whose reads each carry a different intron chain get no transcript at the 2-read
# floor (row 1865: PLA2G4B, 22 reads, 22 chains). This stage (1) assembles again at a 1-read floor ($OUT.multicopy.d/f1.gtf;
# RUSTLE_MULTICOPY_F1=FILE reuses a floor-1 GTF of the SAME bam and seeding pool, e.g. PREFIX.oneread.d/f1.gtf), (2) takes its
# spliced transcripts at places $OUT.gtf leaves empty on that strand, (3) runs the families stage on $OUT.gtf + those candidates
# ($OUT.multicopy.d/u.*) and (4) keeps the candidates whose locus is a copy of a family of >= 2 loci. Genome-wide the candidates
# are 3 to 7 times less often annotated than the rest (row 1867), hence the confinement. $OUT.gtf and every other product are
# untouched; a seeding pool other than `good` is allowed and said.
stage_multicopy() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] multicopy needs $OUT.gtf (run assemble --seed-pool good first), --bam and --fasta" >&2; exit 2; }
  local py=${RUSTLE_PYTHON:-python3} d="$OUT.multicopy.d" f1; mkdir -p "$d"; f1="$d/f1.gtf"
  # RUSTLE_MULTICOPY_MONO=1 (opt-in, row 1919): single-exon models are rescued the same way. The floor-1 run is made without the
  # single-exon support floor (RUSTLE_POLISH_MONO_QUANTILE=0; the shadow and internal-priming rules stay), its single-exon
  # transcripts at places the assembly leaves empty are candidates, and one is kept only as a copy of a family of >= 2 loci
  # (tagged multicopy_rescue "mono"): relatedness to another assembled locus replaces the read-count floor.
  local monoenv=() monoarg=(); case "${RUSTLE_MULTICOPY_MONO:-0}" in 0) ;; 1) monoenv=(RUSTLE_POLISH_MONO_QUANTILE=0); monoarg=(--mono) ;;
    *) echo "[rustle_pipeline] RUSTLE_MULTICOPY_MONO must be 0 or 1 (got '$RUSTLE_MULTICOPY_MONO')" >&2; exit 2 ;; esac
  [ "$SEED_POOL" = good ] || say "multicopy: seeding pool is $SEED_POOL, not good: the floor-1 run will not see near-tie secondaries"
  if [ -n "${RUSTLE_MULTICOPY_F1:-}" ]; then
    [ -s "$RUSTLE_MULTICOPY_F1" ] || { echo "[rustle_pipeline] RUSTLE_MULTICOPY_F1=$RUSTLE_MULTICOPY_F1: no such non-empty file" >&2; exit 2; }
    f1=$RUSTLE_MULTICOPY_F1; say "multicopy: floor-1 assembly taken from $f1"
  else
    # the floor-1 run of stage oneread, with its pooling-guard bypass (RUSTLE_GATE_MIN_READS=1, row 1866). Rows 1865 / 1867 / 1869
    # judged the mode WITHOUT the bypass; row 1937 judged it with: same units and labels on planted reads, more annotated copies at
    # equal junction precision on OR6737 and KB3781, 3-4 x the rescued loci (all labelled).
    # Since 2026-11-24 (user, register 1937) RUSTLE_MULTICOPY_GATE is 1 when unset: the run is then the one-read stage's, and one
    # of the two is reused (f1_shared above). 0 = the run without the bypass, as rows 1865-1869 judged the mode.
    local gate1=(); case "${RUSTLE_MULTICOPY_GATE:-1}" in 0) ;; 1) gate1=(RUSTLE_GATE_MIN_READS=1) ;;
      *) echo "[rustle_pipeline] RUSTLE_MULTICOPY_GATE must be 0 or 1 (got '$RUSTLE_MULTICOPY_GATE')" >&2; exit 2 ;; esac
    local f1set f1have; f1set=$(f1_settings "${RUSTLE_MULTICOPY_GATE:-1}" "${RUSTLE_MULTICOPY_MONO:-0}"); f1have=$(f1_shared "$OUT.oneread.d" "$f1set")
    # (register 1935: giving this run only the alignments at empty places is NOT equivalent, 63% extra candidates; not offered)
    if [ -n "$f1have" ]; then
      f1=$f1have; : > "$OUT.multicopy.log"; say "multicopy: floor-1 assembly taken from $f1 (the one-read stage made it with the same settings)"
    else
      rm -f "$d/f1.gtf" "$d/f1.settings"
      env "${gate1[@]}" "${monoenv[@]}" RUSTLE_DUP_COUNT="${RP_DUP_COUNT_USER:-0}" RUSTLE_MIN_CHAIN_READS=1 RUSTLE_ASSEMBLY_LAYERS=off RUSTLE_ONEREAD=0 RUSTLE_SEED_POOL="$SEED_POOL" RUSTLE_SEED_AS_RATIO="$SEED_RATIO" \
          bash "$0" assemble --bam "$BAM" --fasta "$FASTA" --out "$d/f1" --threads "$THREADS" ${BIN:+--bin "$BIN"} ${CONTIG:+--contig "$CONTIG"} ${AS_TABLE:+--as-table "$AS_TABLE"} > "$OUT.multicopy.log" 2>&1 \
        || { echo "[rustle_pipeline] multicopy: the floor-1 assembly failed (see $OUT.multicopy.log); $OUT.gtf is unaffected" >&2; exit 1; }
      printf '%s\n' "$f1set" > "$d/f1.settings"
    fi
  fi
  $py "$RP_DRIVER_DIR/multicopy_rescue.py" candidates "$OUT.gtf" "$f1" "$d/cand.gtf" "${monoarg[@]}" >> "$OUT.multicopy.log" 2>&1 \
    || { echo "[rustle_pipeline] multicopy: candidates failed (see $OUT.multicopy.log)" >&2; exit 1; }
  say "multicopy: $(tail -1 "$OUT.multicopy.log")"
  rm -f "$d"/u.fam.* "$d/u.families.gtf"
  cat "$OUT.gtf" "$d/cand.gtf" > "$d/u.gtf"
  local ufam="$d/u.gtf"; if [ "$FAM_GTF" != "$OUT.gtf" ]; then cat "$FAM_GTF" "$d/cand.gtf" > "$d/u.families.gtf"; ufam="$d/u.families.gtf"; fi
  ( RUSTLE_FAMREPS=0; RUSTLE_SDCORE=0; FAM_GTF=$ufam; OUT=$d/u; stage_families ) >> "$OUT.multicopy.log" 2>&1   # no `|| exit 1`: it switches set -e off inside
  [ -s "$d/u.fam.copies.tsv" ] || { echo "[rustle_pipeline] multicopy: the families run on the candidates failed (see $OUT.multicopy.log); $OUT.gtf is unaffected" >&2; exit 1; }
  $py "$RP_DRIVER_DIR/multicopy_rescue.py" select "$OUT.gtf" "$d/cand.gtf" "$d/u.fam.copies.tsv" "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.tsv.tmp" >> "$OUT.multicopy.log" 2>&1 \
    || { echo "[rustle_pipeline] multicopy: select failed (see $OUT.multicopy.log)" >&2; rm -f "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.tsv.tmp"; exit 1; }
  mv -f "$OUT.multicopy.tsv.tmp" "$OUT.multicopy.tsv"; mv -f "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.gtf"
  say "multicopy: $(tail -1 "$OUT.multicopy.log") -> $OUT.multicopy.gtf"
  # Counts that do not multiply over indistinguishable copies (row 1871; RUSTLE_MULTICOPY_EQ=0 skips it). A secondary puts one
  # molecule on every near-identical copy, so `reads` is the same molecules at each (three identical TSPY10 copies: 72, 72, 72).
  # eq_classes.py on the mode's GTF (ties on minimap2's ms, exact band) -> $OUT.multicopy.eq.loci.tsv / .classes.tsv; the
  # transcript lines of $OUT.multicopy.gtf get reads_placed / reads_em / reads_max / tie_group / reads_tieadj, and
  # $OUT.multicopy.groups.tsv lists the tie groups with the copies the scores allow (min_expressed .. max_expressed).
  [ "${RUSTLE_MULTICOPY_EQ:-1}" != 0 ] || return 0
  local tab=$OUT.molecules.tsv bed=${RUSTLE_EXCLUDE_BED-$FASTA.exclude.bed}; [ -s "$tab" ] || tab=$BAM.as_table.tsv
  [ -s "$tab" ] || { say "multicopy: no score table ($OUT.molecules.tsv or $BAM.as_table.tsv): tie-adjusted counts skipped"; return 0; }
  local excl=(); [ ! -s "$bed" ] || excl=(--exclude "$bed")
  { $py "$RP_DRIVER_DIR/eq_classes.py" "$OUT.multicopy.gtf" "$BAM" "$tab" "$OUT.multicopy.eq.tmp" --window "${RUSTLE_EQ_WINDOW:-10}" --threads "$JOB_THREADS" "${excl[@]}" \
      && $py "$RP_DRIVER_DIR/multicopy_rescue.py" annotate "$OUT.multicopy.gtf" "$OUT.multicopy.eq.tmp.loci.tsv" "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.groups.tsv.tmp"; } >> "$OUT.multicopy.log" 2>&1 \
    || { echo "[rustle_pipeline] multicopy: tie-adjusted counts failed (see $OUT.multicopy.log); $OUT.multicopy.gtf is kept without them" >&2; rm -f "$OUT.multicopy.eq.tmp".* "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.groups.tsv.tmp"; exit 1; }
  mv -f "$OUT.multicopy.eq.tmp.loci.tsv" "$OUT.multicopy.eq.loci.tsv"; mv -f "$OUT.multicopy.eq.tmp.classes.tsv" "$OUT.multicopy.eq.classes.tsv"
  mv -f "$OUT.multicopy.groups.tsv.tmp" "$OUT.multicopy.groups.tsv"; mv -f "$OUT.multicopy.gtf.tmp" "$OUT.multicopy.gtf"
  say "multicopy: $(tail -1 "$OUT.multicopy.log") -> $OUT.multicopy.groups.tsv"
}

# siblings: OPT-IN, named only (register 1628; never part of `all`). $OUT.sibling.gtf = $OUT.gtf + the chains a family sibling's
# model projects onto a member that exactly one primary alignment carries, from a read with no rival placement (AS table).
# Quality bar passed on three samples (correct share above the 2-read chains), yield bar failed (12-20 correct chains a
# sample): the user ported it as opt-in 2026-11-11. $OUT.gtf and every other product are untouched. Needs python3 + pysam.
stage_siblings() {
  local f; for f in gtf fam.copies.tsv molecules.tsv; do
    [ -s "$OUT.$f" ] || { echo "[rustle_pipeline] siblings needs $OUT.$f (run assemble and families first)" >&2; exit 2; }
  done
  # shortcut: RUSTLE_PYTHON runs UNQUOTED on purpose, so it can be a command with arguments ("pixi run python", the tests' "bash py.sh");
  # a python path that contains a space is therefore not supported (audit 2026-11-19, finding 15 withdrawn)
  local py=${RUSTLE_PYTHON:-python3} b="$RP_DRIVER_DIR/../bench"
  { $py "$b/project/project_models.py" "$OUT.gtf" "$OUT.fam.copies.tsv" "$FASTA" "$BAM" "$OUT.sibling.d" --threads "$THREADS" --as-table "$OUT.molecules.tsv" \
      && $py "$b/famcmp/sibling_chains.py" emit "$OUT.gtf" "$OUT.sibling.d" "$OUT.sibling.gtf.tmp"; } > "$OUT.sibling.log" 2>&1 \
    || { echo "[rustle_pipeline] siblings failed (see $OUT.sibling.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.sibling.gtf.tmp"; exit 1; }
  mv -f "$OUT.sibling.gtf.tmp" "$OUT.sibling.gtf"
  say "siblings: $(tail -1 "$OUT.sibling.log") -> $OUT.sibling.gtf"
}
# flairlike: OPT-IN, named only (register 1640 E5; never part of `all`). $OUT.flairlike.gtf = the chain set FLAIR collapse
# reports (Jaccard 0.87 with FLAIR's spliced chains on gorilla and A119b): every primary alignment's chain as a candidate,
# FLAIR's end choice and subset filter, reads realigned to the candidate sequences, 3-read floor. Costs a second minimap2
# pass (0.5-2.8 h) and precision (half the default's on gorilla, a third on A119b); $OUT.gtf and every other product are untouched.
# Needs python3 + pysam, samtools, minimap2.
stage_flairlike() {
  [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] flairlike needs --bam and --fasta" >&2; exit 2; }
  local py=${RUSTLE_PYTHON:-python3} s="$RP_DRIVER_DIR/../bench/famcmp/flair_like.py" d="$OUT.flairlike.d"
  mkdir -p "$d"
  { FLAIR_ENDS=2 $py "$s" firstpass "$BAM" "$FASTA" "$d/fp" "${MM_THREADS:-$THREADS}" \
      && samtools fastq -F 0x900 -@ 4 "$BAM" | minimap2 -a -t "${MM_THREADS:-$THREADS}" -N 4 "$d/fp.fa" - | $py "$s" gtf "$d/fp.tsv" "$OUT.flairlike.gtf.tmp"; } > "$OUT.flairlike.log" 2>&1 \
    || { echo "[rustle_pipeline] flairlike failed (see $OUT.flairlike.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.flairlike.gtf.tmp"; exit 1; }
  mv -f "$OUT.flairlike.gtf.tmp" "$OUT.flairlike.gtf"
  say "flairlike: $(tail -1 "$OUT.flairlike.log") -> $OUT.flairlike.gtf"
}
# famreps: part of `families` by default since 2026-11-13 (user decision; RUSTLE_FAMREPS=0 = off, the table as before;
# register 1662, rows 1575/1577/1654). Genes with no locus of their own
# (fused or nested in a neighbour's locus) get no family node. Here every exon-disjoint transcript group of a locus with
# >= 2 reads becomes its own representative (locus_reps.py --disjoint), the families stage runs again on that input with
# the same flags (into $OUT.famreps.fam.*), and the extra representatives are attached to the UNCHANGED families, only
# when they cluster beside an unsplit family member (piece_attach.py --keep-all --direct).
# A119b 0.918 / 0.782 -> 0.907 / 0.809 (+58 true, +27 false pairs, 5 of them with no shared exon sequence); HG002
# 0.703 / 0.355 -> 0.727 / 0.444. It FAILS row 1654's A119b bar (precision >= baseline) as scored; default by the user's
# decision. Cost: a second all-vs-all (about the wall of `families`). The first run's table is kept as
# $OUT.fam.base.copies.{tsv,fa,regions}; $OUT.fam.copies.* (what assign reads) holds the extra representatives.
# the extra representatives are built by Python helpers: refuse before the all-vs-all, not after it
famreps_py_guard() {
  [ "${RUSTLE_FAMREPS:-1}" = 0 ] || command -v "${RUSTLE_PYTHON:-python3}" > /dev/null \
    || { echo "[rustle_pipeline] no ${RUSTLE_PYTHON:-python3} on PATH (families adds extra locus representatives with it); set RUSTLE_PYTHON or PATH, or RUSTLE_FAMREPS=0 for the one-representative table" >&2; exit 2; }
}
stage_famreps() {
  famreps_py_guard
  # the table the first families run wrote: PREFIX.fam.base.copies.* once this stage has installed its own as PREFIX.fam.copies.*
  local A="$OUT.fam.copies"; [ -s "$OUT.fam.base.copies.tsv" ] && A="$OUT.fam.base.copies"
  [ -s "$A.tsv" ] && [ -s "$A.fa" ] || { echo "[rustle_pipeline] famreps needs $OUT.fam.copies.tsv/.fa (run families first)" >&2; exit 2; }
  local py=${RUSTLE_PYTHON:-python3} d="$RP_DRIVER_DIR" r="$OUT.famreps" e   # the four helpers sit beside the driver, like annotated_families.py
  $py "$d/locus_reps.py" "$FAM_GTF" "$r.families.gtf" --disjoint > "$OUT.famreps.log" 2>&1 \
    || { echo "[rustle_pipeline] famreps failed (see $OUT.famreps.log); $A.tsv is unaffected" >&2; exit 1; }
  ( RUSTLE_FAMREPS=0; RUSTLE_SDCORE=0; FAM_GTF=$r.families.gtf; OUT=$r; stage_families )  # no `|| exit 1`: it switches set -e off inside
  { $py "$d/piece_merge.py" "$r" "$r.merged.copies.tsv" \
      && $py "$d/piece_attach.py" "$A.tsv" "$r" "$r.merged.copies.tsv" "$r.attached.tsv" --keep-all --direct \
      && $py "$d/reps_table.py" "$A" "$r.fam.copies" "$r.merged.copies.tsv" "$r.attached.tsv" "$r.table" \
      && $py "$d/copy_table_check.py" "$r.table"; } >> "$OUT.famreps.log" 2>&1 \
    || { echo "[rustle_pipeline] famreps failed (see $OUT.famreps.log); $A.tsv is unaffected" >&2; exit 1; }
  if [ "$A" = "$OUT.fam.copies" ]; then for e in tsv fa regions; do mv -f "$OUT.fam.copies.$e" "$OUT.fam.base.copies.$e"; done; fi
  for e in tsv fa regions; do mv -f "$r.table.$e" "$OUT.fam.copies.$e"; done
  say "famreps: $(awk 'NR>1' "$OUT.fam.base.copies.tsv" | wc -l) -> $(awk 'NR>1' "$OUT.fam.copies.tsv" | wc -l) copies ($OUT.fam.copies.tsv; the first run's table is $OUT.fam.base.copies.tsv)"
}
# chainscore: the named stage (register 1761); the same filter is applied to $OUT.final.gtf by default, see final_chainscore. $OUT.chainscore.gtf = $OUT.gtf minus the multi-exon chains the
# frozen logistic score (tools/data/chainscore_*.json; the default model and its cut are named in final_chainscore) puts below the cut; $OUT.gtf is not touched.
# HG002 held out: keeps 0.972 of the correct chains, precision +0.026 over a read-count cut at equal correct chains.
stage_chainscore() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] chainscore needs $OUT.gtf (run assemble first), --bam and --fasta" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/chainscore.py" apply "$(chainscore_model)" "$OUT.gtf" "$BAM" "$FASTA" "$OUT.chainscore.gtf.tmp" "$THREADS" > "$OUT.chainscore.log" 2>&1 \
    || { echo "[rustle_pipeline] chainscore failed (see $OUT.chainscore.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.chainscore.gtf.tmp"; exit 1; }
  mv -f "$OUT.chainscore.gtf.tmp" "$OUT.chainscore.gtf"
  say "chainscore: $(tail -1 "$OUT.chainscore.log") -> $OUT.chainscore.gtf"
}
# trunc5: OPT-IN, named only (register 1677, never part of `all`). $OUT.trunc5.gtf = $OUT.gtf minus 5' truncation chains;
# $OUT.gtf is not touched. Failed its held-out bar on KB3781 (+0.013 over the read-count curve), so it is not a default.
stage_trunc5() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] || { echo "[rustle_pipeline] trunc5 needs $OUT.gtf (run assemble first) and --bam" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/trunc5.py" "$OUT.gtf" "$BAM" "$OUT.trunc5.gtf.tmp" "$THREADS" > "$OUT.trunc5.log" 2>&1 \
    || { echo "[rustle_pipeline] trunc5 failed (see $OUT.trunc5.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.trunc5.gtf.tmp"; exit 1; }
  mv -f "$OUT.trunc5.gtf.tmp" "$OUT.trunc5.gtf"
  say "trunc5: $(tail -1 "$OUT.trunc5.log") -> $OUT.trunc5.gtf"
}
# shortexons: OPT-IN, named only (register 1753 / 1754 / 1757, never part of `all`). minimap2 drops exons of <= 15 bp (row 1746); the
# reads that skip one are rewritten through it (bench/famcmp/short_exons.py) and assembled with the default command into
# $OUT.shortexons.gtf (+ .final.gtf). Measured on eight samples: false skip chains at the reached exons nearly vanish, chains with the
# exon double to triple, panel precision -0.0016 to +0.0049 with no consistent sign (the rewrite also adds new unannotated chains).
# Not a default until a rule removes those. $OUT.gtf is never changed. RUSTLE_SHORTEXONS_KEEP=1 keeps the rewritten BAM.
stage_shortexons() {
  [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] shortexons needs --bam and --fasta" >&2; exit 2; }
  local d="$OUT.shortexons.d"; mkdir -p "$d"
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/short_exons.py" "$BAM" "$FASTA" "$d/reads.bam" "$OUT.shortexons" "$THREADS" > "$OUT.shortexons.log" 2>&1 \
    && RUSTLE_SHORTEXONS=0 bash "$0" assemble --bam "$d/reads.bam" --fasta "$FASTA" --out "$d/run" --threads "$THREADS" --bin "$BIN" > "$OUT.shortexons.assemble.log" 2>&1 \
    && [ -s "$d/run.gtf" ] \
    || { [ -n "${RUSTLE_SHORTEXONS_KEEP:-}" ] || rm -rf "${d:?}"; echo "[rustle_pipeline] shortexons failed (see $OUT.shortexons.log, $OUT.shortexons.assemble.log); $OUT.gtf is unaffected" >&2; exit 1; }
  rm -f "$OUT.shortexons.final.gtf"; mv -f "$d/run.gtf" "$OUT.shortexons.gtf"; [ ! -s "$d/run.final.gtf" ] || mv -f "$d/run.final.gtf" "$OUT.shortexons.final.gtf"
  [ -n "${RUSTLE_SHORTEXONS_KEEP:-}" ] || rm -rf "${d:?}"
  say "shortexons: $(tail -1 "$OUT.shortexons.log") -> $OUT.shortexons.gtf"
}
# capstarts: register 1717 / 1720; `all` writes it last as a side product since 2026-11-14 (user; RUSTLE_CAPSTARTS=0 = off). $OUT.capstarts.gtf = $OUT.gtf plus nested starts
# (>= 5 reads, 3'-terminal sub-run of an emitted chain) with >= 3 capped reads (non-templated 5' G) at one position.
# Refuses libraries without a usable cap G (self-check in the log). $OUT.gtf is never changed.
stage_capstarts() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] capstarts needs $OUT.gtf (run assemble first), --bam and --fasta" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/capg.py" stage "$OUT.gtf" "$BAM" "$FASTA" "$OUT.capstarts.gtf.tmp" "$THREADS" > "$OUT.capstarts.log" 2>&1 \
    || { echo "[rustle_pipeline] capstarts wrote nothing: $(tail -1 "$OUT.capstarts.log"); $OUT.gtf is unaffected" >&2; rm -f "$OUT.capstarts.gtf.tmp"; exit 1; }
  mv -f "$OUT.capstarts.gtf.tmp" "$OUT.capstarts.gtf"
  say "capstarts: $(tail -1 "$OUT.capstarts.log") -> $OUT.capstarts.gtf"
}
# shiftdrop: OPT-IN, named only (register 1728-1730 / 1732 / 1734 / 1738 / 1739, never part of `all`). $OUT.shiftdrop.gtf = $OUT.gtf minus
# weaker shifted-site chains with indels beside the moved intron (bench/famcmp/shiftdrop.py). $OUT.gtf is never changed.
stage_shiftdrop() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] || { echo "[rustle_pipeline] shiftdrop needs $OUT.gtf (run assemble first) and --bam" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/shiftdrop.py" "$OUT.gtf" "$BAM" "$OUT.shiftdrop.gtf.tmp" "$THREADS" > "$OUT.shiftdrop.log" 2>&1 \
    || { echo "[rustle_pipeline] shiftdrop failed (see $OUT.shiftdrop.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.shiftdrop.gtf.tmp"; exit 1; }
  mv -f "$OUT.shiftdrop.gtf.tmp" "$OUT.shiftdrop.gtf"
  say "shiftdrop: $(tail -1 "$OUT.shiftdrop.log") -> $OUT.shiftdrop.gtf"
}
# earlyends: register 1721 / 1722; `all` writes it last as a side product since 2026-11-14 (RUSTLE_EARLYENDS=0 = off). $OUT.earlyends.gtf = $OUT.gtf plus nested early ends
# (>= 5 reads, 5'-terminal sub-run of an emitted chain) whose reads end at one position with a polyA signal and no A-rich run.
stage_earlyends() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] earlyends needs $OUT.gtf (run assemble first), --bam and --fasta" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/capg.py" stage "$OUT.gtf" "$BAM" "$FASTA" "$OUT.earlyends.gtf.tmp" "$THREADS" --ends > "$OUT.earlyends.log" 2>&1 \
    || { echo "[rustle_pipeline] earlyends failed (see $OUT.earlyends.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.earlyends.gtf.tmp"; exit 1; }
  mv -f "$OUT.earlyends.gtf.tmp" "$OUT.earlyends.gtf"
  say "earlyends: $(tail -1 "$OUT.earlyends.log") -> $OUT.earlyends.gtf"
}
# relatives: OPT-IN, named only (register 1676 / 1737, never part of `all`). Reads the families stage output; changes none of it.
stage_relatives() {
  local f s="$RP_DRIVER_DIR/../bench/denovo_eval/relative_links.py" bed=${RUSTLE_SEDEF_BED:-$FASTA.sedef.bed}
  for f in loci.gff3 loci.gated.paf copies.tsv; do
    [ -s "$OUT.fam.$f" ] || { echo "[rustle_pipeline] relatives needs $OUT.fam.$f (run families first)" >&2; exit 2; }
  done
  [ -s "$bed" ] || { echo "[rustle_pipeline] relatives needs a SEDEF BED of this genome (RUSTLE_SEDEF_BED or $FASTA.sedef.bed): SD links are part of the score" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$s" --sedef "$bed" ${RUSTLE_EXCLUDE_BED:+--exclude-bed "$RUSTLE_EXCLUDE_BED"} "$OUT.fam" "$OUT.fam.relatives.tsv" > "$OUT.fam.relatives.log" 2>&1 \
    || { echo "[rustle_pipeline] relatives failed (see $OUT.fam.relatives.log); the families output is unaffected" >&2; exit 1; }
  say "relatives: $(tail -1 "$OUT.fam.relatives.log") (table $OUT.fam.relatives.tsv; links, not members)"
}
# sdcore: register 1755 / 1756; `families` runs it last when a SEDEF BED is set (default since 2026-11-15, user; RUSTLE_SDCORE=0 = off).
# Reads PREFIX.fam.copies.tsv; changes nothing.
stage_sdcore() {
  local s="$RP_DRIVER_DIR/../bench/denovo_eval/sdcore.py" bed=${RUSTLE_SEDEF_BED:-$FASTA.sedef.bed}
  [ -s "$OUT.fam.copies.tsv" ] || { echo "[rustle_pipeline] sdcore needs $OUT.fam.copies.tsv (run families first)" >&2; exit 2; }
  [ -s "$bed" ] || { echo "[rustle_pipeline] sdcore needs a SEDEF BED of this genome (RUSTLE_SEDEF_BED or $FASTA.sedef.bed)" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$s" "$OUT.fam.copies.tsv" "$bed" "$OUT.fam.sdcore.tsv" > "$OUT.fam.sdcore.log" 2>&1 \
    || { echo "[rustle_pipeline] sdcore failed (see $OUT.fam.sdcore.log); the families output is unaffected" >&2; exit 1; }
  say "sdcore: $(tail -1 "$OUT.fam.sdcore.log") (table $OUT.fam.sdcore.tsv; a label, the families are unchanged)"
}
# exonmap: OPT-IN, named only (register 1595, never part of `all`). A view of the families stage output; changes none of it.
stage_exonmap() {
  local f s="$RP_DRIVER_DIR/../bench/denovo_eval/exonmap.py"
  for f in clusters.tsv loci.gff3 loci.gated.paf; do
    [ -s "$OUT.fam.$f" ] || { echo "[rustle_pipeline] exonmap needs $OUT.fam.$f (run families first)" >&2; exit 2; }
  done
  rm -rf "$OUT.fam.exonmap.gfa"
  { ${RUSTLE_PYTHON:-python3} "$s" --threads "$THREADS" "$OUT.fam" "$OUT.fam.exonmap" \
      && ${RUSTLE_PYTHON:-python3} "$s" gfa "$OUT.fam" "$OUT.fam.exonmap.columns.tsv" "$OUT.fam.exonmap.gfa"; } > "$OUT.fam.exonmap.log" 2>&1 \
    || { echo "[rustle_pipeline] exonmap failed (see $OUT.fam.exonmap.log); the families output is unaffected" >&2; exit 1; }
  say "exonmap: $(tail -1 "$OUT.fam.exonmap.log") (table $OUT.fam.exonmap.columns.tsv)"
}
# repeatloci: OPT-IN, named only (registers 1698 / 1702, never part of `all`). $OUT.repeat_loci.tsv = one row per assembly
# locus with the soft-masked share of its exons, the share covered by the rDNA unit (RUSTLE_REPEAT_FASTA, default
# tools/data/rdna_human.fa) and repeat = 1 when the first is >= 0.9 or the second >= 0.5. A table only: nothing is dropped.
stage_repeatloci() {
  local ref="${RUSTLE_REPEAT_FASTA:-$RP_DRIVER_DIR/data/rdna_human.fa}"
  [ -s "$OUT.gtf" ] && [ -s "$FASTA" ] && [ -s "$ref" ] || { echo "[rustle_pipeline] repeatloci needs $OUT.gtf (run assemble first), --fasta and $ref" >&2; exit 2; }
  RUSTLE_SEQFAM_THREADS="$THREADS" ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/repeat_loci.py" "$OUT.gtf" "$FASTA" - "$ref" "$OUT.repeat_loci.tmp" > "$OUT.repeat_loci.log" 2>&1 \
    || { echo "[rustle_pipeline] repeatloci failed (see $OUT.repeat_loci.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.repeat_loci.tmp.tsv" "$OUT.repeat_loci.tmp.fa"; exit 1; }
  mv -f "$OUT.repeat_loci.tmp.tsv" "$OUT.repeat_loci.tsv"; rm -f "$OUT.repeat_loci.tmp.fa"
  say "repeatloci: $(tail -1 "$OUT.repeat_loci.log") -> $OUT.repeat_loci.tsv"
}
# altends: OPT-IN, named only (never part of `all`). $OUT.altends.gtf = $OUT.gtf with attributes on multi-exon transcript lines
# whose commonest full-chain read 3' end is not the emitted (most extended) one: alt_tes, alt_tes_reads, alt_tes_polya_signal,
# alt_tes_a_rich. No coordinate changes; $OUT.gtf is untouched. Why: about a quarter of A119b's correct chains mostly use a
# polyA site upstream of the annotated one (docs/STEP0_tool_match_2026-11-11.md §11). Needs python3 + pysam.
stage_altends() {
  [ -s "$OUT.gtf" ] && [ -s "$BAM" ] && [ -s "$FASTA" ] || { echo "[rustle_pipeline] altends needs $OUT.gtf (run assemble first), --bam and --fasta" >&2; exit 2; }
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/../bench/famcmp/alt_ends.py" "$OUT.gtf" "$BAM" "$FASTA" "$OUT.altends.gtf.tmp" "$THREADS" > "$OUT.altends.log" 2>&1 \
    || { echo "[rustle_pipeline] altends failed (see $OUT.altends.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.altends.gtf.tmp"; exit 1; }
  mv -f "$OUT.altends.gtf.tmp" "$OUT.altends.gtf"
  say "altends: $(tail -1 "$OUT.altends.log") -> $OUT.altends.gtf"
}
# psi: OPT-IN, named only (2026-11-16). One BAM pass; $OUT.psi.tsv = per-intron PSI within bunches with tie ranges; nothing reads it.
stage_psi() {
  local tab=$OUT.molecules.tsv
  [ -s "$OUT.gtf" ] || { echo "[rustle_pipeline] psi needs $OUT.gtf (run assemble first)" >&2; exit 2; }
  [ -s "$tab" ] || tab=$BAM.as_table.tsv
  [ -s "$tab" ] || { echo "[rustle_pipeline] psi needs the score table $OUT.molecules.tsv or $BAM.as_table.tsv (assemble writes it)" >&2; exit 2; }
  # with exon columns, also the raw count of every intron: the family view's rows for copies under the 3-read floor
  local famview="" dump=() extra=()
  if [ -s "$OUT.fam.exonmap.columns.tsv" ] && [ -s "$OUT.fam.loci.gff3" ] && [ ! "$OUT.fam.exonmap.columns.tsv" -ot "$OUT.fam.loci.gff3" ]; then
    famview=1; dump=(--introns "$OUT.psi.introns.tsv")
  fi
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/psi_bunches.py" "$OUT.gtf" "$BAM" "$tab" "$OUT.psi.tsv" --window "${RUSTLE_EQ_WINDOW:-10}" ${RUSTLE_EQ_RATIO:+--ratio "$RUSTLE_EQ_RATIO"} --threads "$JOB_THREADS" ${dump[@]+"${dump[@]}"} > "$OUT.psi.log" 2>&1 \
    || { echo "[rustle_pipeline] psi failed (see $OUT.psi.log); no other product is affected" >&2; rm -f "$OUT.psi.tsv.tmp"; exit 1; }
  say "psi: $(tr '\n' ' ' < "$OUT.psi.log")-> $OUT.psi.tsv"
  # the same junction across the copies of a family, when stage `exonmap` has written its columns (else skipped)
  if [ -n "$famview" ]; then
    if [ -s "$OUT.fam.loci.gated.paf" ]; then extra=(--introns "$OUT.psi.introns.tsv" --paf "$OUT.fam.loci.gated.paf")
    else say "psi: no $OUT.fam.loci.gated.paf (the families stage writes it): no added rows for copies under the read floor"; fi
    ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/psi_family.py" "$OUT.psi.tsv" "$OUT.fam.exonmap.columns.tsv" "$OUT.fam.loci.gff3" "$OUT.fam.psi.tsv" ${extra[@]+"${extra[@]}"} > "$OUT.fam.psi.log" 2>&1 \
      || { echo "[rustle_pipeline] psi: the family view failed (see $OUT.fam.psi.log); $OUT.psi.tsv is complete" >&2; rm -f "$OUT.fam.psi.tsv.tmp"; exit 1; }
    say "psi: $(tr '\n' ' ' < "$OUT.fam.psi.log")-> $OUT.fam.psi.tsv"
  else
    rm -f "$OUT.fam.psi.tsv"
    say "psi: no $OUT.fam.exonmap.columns.tsv newer than $OUT.fam.loci.gff3 (run stage exonmap), family view skipped"
  fi
}
# eqclasses: OPT-IN, named only (2026-11-16). One BAM pass over the loci `families` reads; nothing reads its tables.
# Ties are read on minimap2's ms (not AS, which ties spliced genes with processed pseudogenes); the score table needs its ms columns.
stage_eqclasses() {
  fam_gtf_guard eqclasses
  local tab=$OUT.molecules.tsv gtf=$FAM_GTF bed=${RUSTLE_EXCLUDE_BED-$FASTA.exclude.bed}   # set empty = off, as in assign
  local excl=(); [ ! -s "$bed" ] || excl=(--exclude "$bed")
  local fw=(); [ -z "${RP_FAMWISE:-}" ] || fw=(--famwise "$OUT.eq.tmp.fw")   # stage famwise: the same pass also writes its two tables
  [ -s "$tab" ] || tab=$BAM.as_table.tsv
  [ -s "$tab" ] || { echo "[rustle_pipeline] eqclasses needs the score table $OUT.molecules.tsv or $BAM.as_table.tsv (assemble writes it)" >&2; exit 2; }
  if [ -n "${NAMES:-}" ]; then
    ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/gtf_names.py" "$NAMES" "$FAM_GTF" "$OUT.eq.named.gtf" > "$OUT.eq.log" 2>&1 \
      || { echo "[rustle_pipeline] eqclasses: naming failed (see $OUT.eq.log)" >&2; exit 1; }
    gtf=$OUT.eq.named.gtf
  fi
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/eq_classes.py" "$gtf" "$BAM" "$tab" "$OUT.eq.tmp" --window "${RUSTLE_EQ_WINDOW:-10}" ${RUSTLE_EQ_RATIO:+--ratio "$RUSTLE_EQ_RATIO"} --threads "$JOB_THREADS" \
    "${excl[@]}" "${fw[@]}" >> "$OUT.eq.log" 2>&1 \
    || { echo "[rustle_pipeline] eqclasses failed (see $OUT.eq.log); no other product is affected" >&2; rm -f "$OUT.eq.tmp".*; exit 1; }
  mv -f "$OUT.eq.tmp.loci.tsv" "$OUT.eq.loci.tsv"; mv -f "$OUT.eq.tmp.classes.tsv" "$OUT.eq.classes.tsv"
  say "eqclasses: $(tail -2 "$OUT.eq.log" | tr '\n' ' ')-> $OUT.eq.loci.tsv, $OUT.eq.classes.tsv"
}
# famwise: OPT-IN, named only (2026-11-18; docs/superpowers/specs/2026-11-18-family-wise-assembly-design.md). `eqclasses` with its
# --famwise switch (one BAM pass; the eq tables are byte-identical), then tools/famwise.py on the tables. Nothing reads its products.
stage_famwise() {
  [ -z "${RUSTLE_EQ_RATIO:-}" ] || { echo "[rustle_pipeline] famwise needs the window band: unset RUSTLE_EQ_RATIO (RUSTLE_EQ_WINDOW sets the band)" >&2; exit 2; }
  RP_FAMWISE=1 stage_eqclasses
  local gtf=$FAM_GTF fw=$OUT.famwise
  [ -z "${NAMES:-}" ] || gtf=$OUT.eq.named.gtf
  # the two tables of the BAM pass are kept whatever the layer does; its own four are never left stale or partial
  mv -f "$OUT.eq.tmp.fw.reads.tsv" "$fw.reads.tsv"; mv -f "$OUT.eq.tmp.fw.placed.tsv" "$fw.placed.tsv"
  rm -f "$fw.tx.tsv" "$fw.classes.tsv" "$fw.families.tsv" "$fw.gtf"
  ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/famwise.py" "$gtf" "$fw" --eq-loci "$OUT.eq.loci.tsv" --window "${RUSTLE_EQ_WINDOW:-10}" --min-link "${RUSTLE_FW_MIN_LINK:-1}" > "$OUT.famwise.log" 2>&1 \
    || { echo "[rustle_pipeline] famwise failed (see $OUT.famwise.log); $fw.reads.tsv and .placed.tsv are kept, no other product is affected" >&2; rm -f "$fw.tx.tsv" "$fw.classes.tsv" "$fw.families.tsv" "$fw.gtf"; exit 1; }
  say "famwise: $(tail -3 "$OUT.famwise.log" | tr '\n' ' ')-> $OUT.famwise.families.tsv, .tx.tsv, .classes.tsv, .reads.tsv, .gtf"
}
# names: --names ANNOT (opt-in). Copies of the assemblies with the annotation's names as extra attributes; nothing reads them.
stage_names() {
  [ -s "$NAMES" ] || { echo "[rustle_pipeline] names needs --names ANNOT.gff3|gtf" >&2; exit 2; }
  [ -s "$OUT.gtf" ] || { echo "[rustle_pipeline] names needs $OUT.gtf (run assemble first)" >&2; exit 2; }
  local g
  for g in "$OUT" "$OUT.final"; do
    [ -s "$g.gtf" ] || continue
    ${RUSTLE_PYTHON:-python3} "$RP_DRIVER_DIR/gtf_names.py" "$NAMES" "$g.gtf" "$g.named.gtf.tmp" > "$OUT.names.log" 2>&1 \
      || { echo "[rustle_pipeline] names failed (see $OUT.names.log); $g.gtf is unaffected" >&2; rm -f "$g.named.gtf.tmp"; exit 1; }
    mv -f "$g.named.gtf.tmp" "$g.named.gtf"
    say "names: $(tail -1 "$OUT.names.log") -> $g.named.gtf"
  done
}
stage_guided() {
  [ -n "$GUIDE" ] || { echo "[rustle_pipeline] guided needs --guide GTF" >&2; exit 2; }
  [ -s "$OUT.gtf" ] || { echo "[rustle_pipeline] guided needs $OUT.gtf (run assemble first)" >&2; exit 2; }
  # register 1934 (user 2026-11-24): the rule reads the layered, chain-score-filtered assembly $OUT.final.gtf (10 of 10 judging
  # samples: fewer wrong chains, 1-31 hidden transcripts lost); $OUT.gtf when there is no layered assembly, or RUSTLE_GUIDED_INPUT=raw
  local gin="$OUT.final.gtf"
  if [ "${RUSTLE_GUIDED_INPUT:-final}" = raw ] || [ ! -s "$gin" ]; then gin="$OUT.gtf"; fi
  case "${RUSTLE_GUIDED_INPUT:-final}" in final|raw) ;; *) echo "[rustle_pipeline] RUSTLE_GUIDED_INPUT must be final or raw" >&2; exit 2;; esac
  say "guided: rule applied to $gin"
  "$BIN/utilities" guide --gtf "$gin" --guide "$GUIDE" --bam "$BAM" --threads "$THREADS" --out "$OUT.guided.gtf.tmp" 2> "$OUT.guided.log" \
    || { echo "[rustle_pipeline] --guide failed (see $OUT.guided.log); $OUT.gtf is unaffected" >&2; rm -f "$OUT.guided.gtf.tmp"; exit 1; }
  mv -f "$OUT.guided.gtf.tmp" "$OUT.guided.gtf"
  say "guided: $(awk -F'\t' '$3=="transcript"' "$OUT.guided.gtf" | wc -l) transcripts in the guided assembly $OUT.guided.gtf ($(tail -1 "$OUT.guided.log" | cut -c1-300))"
}
# families: the de novo families AND their copy table (--emit-units with --from-gtf: PREFIX.fam.copies.tsv/.fa/.regions,
# the gw_family_catalog copies contract, one copy per member locus = its representative transcript and its spliced exon
# sum). clusters.tsv / loci.* are byte-identical to a run without it (cmp-checked 2026-09-25,
# docs/archive/2026-09/PREREG_families_copy_table_2026-09-25.md); the copy table is a new product. A binary older than the copy table
# (its --help does not name <out>.copies.tsv) still writes the families, with a warning and no copy table.
stage_families() {
  fam_gtf_guard families
  famreps_py_guard
  # RUSTLE_FAMILIES_EXCLUDE_BED=BED (OPT-IN, register 1749; unset = the same command as before): loci whose exon bases lie mostly
  # inside an element of the BED (T2T censat satellites / rDNA) are left out of the families input, so they never become family
  # nodes (assign's RUSTLE_EXCLUDE_BED skips such families after they are built). The filtered input is PREFIX.fam.input.gtf;
  # the dropped loci are in PREFIX.fam.input.gtf.excluded.txt. `local`: famreps, called below, reads the same filtered input.
  if [ -n "${RUSTLE_FAMILIES_EXCLUDE_BED:-}" ]; then
    [ -s "$RUSTLE_FAMILIES_EXCLUDE_BED" ] || { echo "[rustle_pipeline] RUSTLE_FAMILIES_EXCLUDE_BED=$RUSTLE_FAMILIES_EXCLUDE_BED: no such non-empty file" >&2; exit 2; }
    "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/gtf_exclude_bed.py" "$FAM_GTF" "$RUSTLE_FAMILIES_EXCLUDE_BED" "$OUT.fam.input.gtf" 2> "$OUT.fam.input.log" \
      || { echo "[rustle_pipeline] families: gtf_exclude_bed.py failed (see $OUT.fam.input.log)" >&2; exit 1; }
    say "families: exclude-bed: $(cat "$OUT.fam.input.log")"
    local FAM_GTF=$OUT.fam.input.gtf
  fi
  say "families: gene families on the de novo loci of $FAM_GTF"
  local copies=() help
  help=$("$BIN/mcl_families" --help 2>&1 || true)
  case "$help" in
    *'<out>.copies.tsv'*) copies=(--emit-units);;
    *) say "families: WARNING $BIN/mcl_families predates the families copy table; rebuild it to write $OUT.fam.copies.tsv";;
  esac
  if [ ${#FAM_EXTRA[@]} -gt 0 ] && [[ "$help" != *'--emit-container'* ]]; then
    echo "[rustle_pipeline] RUSTLE_FAMILY_CONTAINER=1 but $BIN/mcl_families predates --emit-container; rebuild it" >&2; exit 2
  fi
  if [ ${#FAM_COV[@]} -gt 0 ] && [[ "$help" != *'--min-cov-shorter'* ]]; then
    echo "[rustle_pipeline] RUSTLE_MIN_COV_SHORTER is set but $BIN/mcl_families predates --min-cov-shorter; rebuild it" >&2; exit 2
  fi
  if [ ${#FAM_REL[@]} -gt 0 ] && [[ "$help" != *'--emit-relations'* ]]; then
    echo "[rustle_pipeline] RUSTLE_FAMILY_RELATIONS=1 but $BIN/mcl_families predates --emit-relations; rebuild it" >&2; exit 2
  fi
  # RUSTLE_REPRESENTATIVE=most-reads|most-junctions (unset = nothing passed: the binary's default, most-reads, the earlier
  # products byte for byte; `most-reads` names that default): mcl_families --representative, which transcript represents each
  # de novo locus (docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md). Read here, in this stage's code alone, so that
  # figures/samples.py's driver_stage_code leaves the other stages' code hashes as they were.
  local rep=()
  case "${RUSTLE_REPRESENTATIVE:-}" in
    "") ;;
    most-reads|most-junctions)
      [[ "$help" == *'--representative'* ]] || { echo "[rustle_pipeline] RUSTLE_REPRESENTATIVE is set but $BIN/mcl_families predates --representative; rebuild it" >&2; exit 2; }
      rep=(--representative "$RUSTLE_REPRESENTATIVE") ;;
    *) echo "[rustle_pipeline] RUSTLE_REPRESENTATIVE must be most-reads or most-junctions (got '$RUSTLE_REPRESENTATIVE')" >&2; exit 2 ;;
  esac
  # Edge rule. DEFAULT since 2026-11-09 (user decision; row 1593 PASS on HG002 held-out): `share` = an edge needs one record on
  # which >= half of the fewer-exon locus's exons match an exon of the other (`--min-exon-share`, RUSTLE_MIN_EXON_SHARE,
  # default 0.5), with the coverage gate off (--min-shared-exon-frac 0.0) and MCL inflation RUSTLE_INFLATION (default 2.0).
  # HG002: plain .545/.396 vs the coverage rule's .466/.345; with RUSTLE_SEDEF_BED .703/.355 vs .719/.321.
  # RUSTLE_FAMILY_EDGE_RULE=coverage = the rule before 2026-11-09: RUSTLE_MIN_SHARED_EXON_FRAC=F (r1528; default 0.60) and the
  # binary's inflation (2.8); setting RUSTLE_MIN_SHARED_EXON_FRAC alone also selects it.
  local edge=()
  case "${RUSTLE_FAMILY_EDGE_RULE:-$([ -n "${RUSTLE_MIN_SHARED_EXON_FRAC:-}" ] && echo coverage || echo share)}" in
    share)
      [[ "$help" == *'--min-exon-share'* ]] || { echo "[rustle_pipeline] the default edge rule needs an mcl_families with --min-exon-share; rebuild it (or RUSTLE_FAMILY_EDGE_RULE=coverage)" >&2; exit 2; }
      edge=(--min-shared-exon-frac 0.0 --min-exon-share "${RUSTLE_MIN_EXON_SHARE:-0.5}" --inflation "${RUSTLE_INFLATION:-2.0}") ;;
    coverage)
      case "${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}" in
        0.[0-9]*|1|1.0) ;;
        *) echo "[rustle_pipeline] RUSTLE_MIN_SHARED_EXON_FRAC must be a fraction 0-1 (got '$RUSTLE_MIN_SHARED_EXON_FRAC')" >&2; exit 2 ;;
      esac
      edge=(--min-shared-exon-frac "${RUSTLE_MIN_SHARED_EXON_FRAC:-0.60}")
      [ -n "${RUSTLE_INFLATION:-}" ] && edge+=(--inflation "$RUSTLE_INFLATION") ;;
    *) echo "[rustle_pipeline] RUSTLE_FAMILY_EDGE_RULE must be share or coverage (got '$RUSTLE_FAMILY_EDGE_RULE')" >&2; exit 2 ;;
  esac
  say "families: edge rule ${edge[*]}"
  # RUSTLE_SEDEF_BED=FILE (SEDEF pairs BED of THIS genome; human CHM13 and gorilla have one, docs/DATA.md): SD damping ON at
  # RUSTLE_SD_DAMP (default 0.25, row 1580 PASS on HG002 held-out, default per user 2026-11-09; 1 = off). Unset = no damping,
  # the same command as before. The factor is always passed, so an older binary cannot fall back to its own default.
  local sd=()
  # the genome's own SEDEF BED as a sidecar, FASTA.sedef.bed, is used when RUSTLE_SEDEF_BED is unset (round-4 audit: SD damping
  # is the default for genomes that have SEDEF, and an unset variable used to drop it silently); RUSTLE_SEDEF_BED= (empty) turns it off
  if [ -z "${RUSTLE_SEDEF_BED+x}" ] && [ -s "$FASTA.sedef.bed" ]; then RUSTLE_SEDEF_BED=$FASTA.sedef.bed; fi
  if [ -n "${RUSTLE_SEDEF_BED:-}" ]; then
    [ -s "$RUSTLE_SEDEF_BED" ] || { echo "[rustle_pipeline] RUSTLE_SEDEF_BED=$RUSTLE_SEDEF_BED: no such non-empty file" >&2; exit 2; }
    [[ "$help" == *'--sd-damp-bed'* ]] || { echo "[rustle_pipeline] RUSTLE_SEDEF_BED is set but $BIN/mcl_families predates --sd-damp-bed; rebuild it" >&2; exit 2; }
    sd=(--sd-damp-bed "$RUSTLE_SEDEF_BED" --sd-damp "${RUSTLE_SD_DAMP:-0.25}")
    say "families: SD damping ${RUSTLE_SD_DAMP:-0.25} with $RUSTLE_SEDEF_BED"
  fi
  # an earlier run's copy table must not survive a run that writes none (a binary older than the copy table): assign would read it
  rm -f "$OUT".fam.copies.{tsv,fa,regions,merged.tsv} "$OUT".fam.base.copies.{tsv,fa,regions} "$OUT.fam.sdcore.tsv"
  # threads: JOB_THREADS (the all-vs-all is ~98% of this stage; A119b 72 threads 410 s, the same families as 12 threads, 2026-11-09)
  local fth=$MM_THREADS
  say "families: $fth threads"
  "$BIN/mcl_families" --from-gtf "$FAM_GTF" --fasta "$FASTA" --threads "$fth" \
    --min-exonic-bp 1 "${edge[@]}" "${copies[@]}" "${FAM_EXTRA[@]}" "${FAM_REL[@]}" "${FAM_COV[@]}" "${rep[@]}" "${sd[@]}" --out "$OUT.fam" > "$OUT.families.log" 2>&1
  # without minimap2-gate on PATH (or RUSTLE_MINIMAP2) every locus pair is aligned: same families, ~5x the wall (audit round 1)
  if grep -q 'has no --gate' "$OUT.families.log"; then
    say "families: WARNING no minimap2-gate on PATH: every pair was aligned (~5x slower); build it with tools/minimap2_gate/build.sh"
  fi
  say "families: $(awk 'NR>1' "$OUT.fam.clusters.tsv" | cut -f1 | sort -u | wc -l) clusters ($OUT.fam.clusters.tsv)"
  if [ ${#copies[@]} -gt 0 ]; then
    say "families: $(awk 'NR>1' "$OUT.fam.copies.tsv" | wc -l) copies (locus representatives) in $(awk 'NR>1' "$OUT.fam.copies.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.fam.copies.tsv)"
  fi
  if [ ${#FAM_REL[@]} -gt 0 ]; then
    say "families: relations $(awk 'NR>1' "$OUT.fam.relations.tsv" | wc -l) split transcripts, $(awk -F'\t' 'NR>1 && $14=="cover"' "$OUT.fam.relations.tsv" | wc -l) cover, $(awk 'NR>1' "$OUT.fam.members_by_locus.tsv" | wc -l) members by locus ($OUT.fam.relations.tsv, $OUT.fam.members_by_locus.tsv)"
  fi
  if [ ${#FAM_EXTRA[@]} -gt 0 ]; then
    say "families: container $(awk -F'\t' '$1=="blocks"{b=$2} $1=="accessory_blocks"{a=$2} $1=="family_relations_directed"{r=$2} END{print b" blocks, "a" accessory, "r" directed family relations"}' "$OUT.fam.container_summary.tsv") ($OUT.fam.container.tsv)"
  fi
  # extra representatives for genes with no locus of their own (stage_famreps; default since 2026-11-13, RUSTLE_FAMREPS=0 = off)
  if [ "${RUSTLE_FAMREPS:-1}" != 0 ] && [ -s "$OUT.fam.copies.tsv" ] && [ -s "$OUT.fam.copies.fa" ]; then stage_famreps; fi
  # sd_core label on the final copy table (stage_sdcore; default since 2026-11-15, user, whenever the genome has a SEDEF BED;
  # RUSTLE_SDCORE=0 = off). A side table: a failure is reported and leaves the families output complete.
  # Row 1756: flagged pairs are enriched for true pairs, NOT confirmed (HG002 0.795 true against its 0.90 bar).
  if [ "${RUSTLE_SDCORE:-1}" != 0 ] && [ -n "${RUSTLE_SEDEF_BED:-}" ] && [ -s "$OUT.fam.copies.tsv" ] && command -v "${RUSTLE_PYTHON:-python3}" > /dev/null; then
    ( stage_sdcore ) || say "families: sdcore wrote no $OUT.fam.sdcore.tsv (see $OUT.fam.sdcore.log); every other product is complete"
  fi
}
# catalog: LEGACY (2026-09-25). gw_family_catalog's copy catalog (primaries only, span-aware POA collapse, exon-sum k11
# edges, gamma-quasi-clique) was the O2 roster before the families stage wrote its own copy table; kept runnable for
# comparison and for the runs made with it. The default de novo family definition is the families stage.
stage_catalog() {
  say "catalog: gw_family_catalog on $BAM"
  if [ "$PIECEWISE" = 1 ]; then
    [ "$CACHE" = 1 ] || { echo "catalog --piecewise keeps its pieces in PREFIX.cache: drop --no-cache" >&2; exit 2; }
    local rc=0
    echo "=== [rustle_pipeline] $(date '+%F %T') catalog --piecewise --max-pieces $MAX_PIECES --budget-s $BUDGET_S --piece-records $PIECE_RECORDS" >> "$OUT.catalog.log"
    env "${INSPECT_CAT[@]}" "$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" --threads "$THREADS" --out "$OUT.cat" \
      --piecewise --max-pieces "$MAX_PIECES" --budget-s "$BUDGET_S" --piece-records "$PIECE_RECORDS" ${PIECE:+--piece "$PIECE"} \
      >> "$OUT.catalog.log" 2>&1 || rc=$?
    if [ "$rc" = 75 ]; then
      say "catalog: $(grep -oE '[0-9]+ of [0-9]+ pieces cached[^;]*|all-vs-all aligner stopped with its progress kept' "$OUT.catalog.log" | tail -1); call again to continue (exit 75)"
      exit 75
    fi
    [ "$rc" = 0 ] || { say "catalog failed (exit $rc), see $OUT.catalog.log"; exit "$rc"; }
  else
    env "${INSPECT_CAT[@]}" "$BIN/gw_family_catalog" --bam "$BAM" --fasta "$FASTA" --threads "$THREADS" --out "$OUT.cat" > "$OUT.catalog.log" 2>&1
  fi
  say "catalog: $(awk 'NR>1' "$OUT.cat.copies.tsv" | wc -l) copies in $(awk 'NR>1 && $2>=2' "$OUT.cat.families.tsv" | wc -l) multi-copy families"
}
# candidates (O3; spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md §4, §7; OPT-IN since 2026-10-02, ruling R14:
# its pre-registered acceptance, Amendment 12, FAILED, docs/archive/2026-10/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md; the re-run, Amendment 13,
# passed, but its no-deletion control, Amendment 14, failed, so the 2026-10-03 default flip was reverted,
# docs/archive/2026-10/O3_CANDIDATES_CONTROL_A14_2026-10-03.md — naming the stage runs it, `all` runs it only with --candidates, and assign
# and flag use its products only with --candidates): candidate_copies turns each family's read net (reads with a record on its
# copies, plus the unmapped and the poorly placed reads >= 300 bp that align, map-ont, over >= 50% of their length at
# de <= 0.20 to the run's net reads or copies: Amendment 13b) into read clusters at --delta, one consensus per cluster on its
# structurally central member (Amendments 13d/13e), and candidate copies (clusters beyond delta of every reference locus,
# merged by the significance test, flagged with >= 6 reads), each represented by the exon union of its clusters ->
# PREFIX.cand.{candidates.tsv,contigs.fa,nets.fa,...}. With a flagged candidate: `utilities candidate-augment` writes PREFIX.aug.{fa,copies.tsv,copies.fa,regions.txt,
# families.txt} (each candidate a contig `cand_<family>_<k>` of PREFIX.aug.fa and a `member_status candidate` row of its
# family), and the reads of the candidate families (PREFIX.cand.nets.fa, every read of each such family's net) are
# realigned to PREFIX.aug.fa with the pipeline's own minimap2 flags -> PREFIX.aug.bam, the reads `assign --candidates` gives
# those families. That realignment is the placement of the candidates' reads (ruling R13): a read whose best alignment is on
# a candidate contig lands there (the fixture: all 60 reads of the deleted copy, MAPQ 60). An earlier run's PREFIX.aug.* are
# removed first, so a run without a flagged candidate leaves none behind. PREFIX.aug.fa is a whole-genome copy that minimap2
# indexes again on every run (about 400 s and 19 GB on a human genome).
stage_candidates() {
  [ -n "$INDEX" ] || { echo "candidates needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  [ -s "$OUT.fam.copies.tsv" ] || { echo "candidates needs $OUT.fam.copies.tsv (run the families stage)" >&2; exit 2; }
  rm -f "$OUT".aug.{fa,fa.fai,copies.tsv,copies.fa,regions.txt,families.txt,bam,bam.bai,mm2.log} "$OUT".cand.candidates.tsv
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    say "candidates: $OUT.fam.copies.tsv lists no copy (no family): nothing to do"; return 0
  fi
  say "candidates: candidate_copies on $OUT.fam.copies.tsv (delta $DELTA, at most $CAND_MAX reads per family)"
  "$BIN/candidate_copies" --bam "$BAM" --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --index "$INDEX" --delta "$DELTA" --max-reads "$CAND_MAX" --threads "$THREADS" --out "$OUT.cand" > "$OUT.candidates.log" 2>&1 \
    || { local rc=$?; say "candidates: candidate_copies failed (exit $rc), see $OUT.candidates.log"; exit "$rc"; }
  local n
  n=$(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | wc -l)
  say "candidates: $n flagged candidate copies in $(awk -F'\t' 'NR>1 && $5==1' "$OUT.cand.candidates.tsv" | cut -f1 | sort -u | wc -l) families ($OUT.cand.candidates.tsv)"
  [ "$n" -gt 0 ] || return 0
  "$BIN/utilities" candidate-augment --fasta "$FASTA" --copies "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --regions "$OUT.fam.copies.regions" --cand "$OUT.cand" --out "$OUT.aug" || exit $?
  samtools faidx "$OUT.aug.fa"
  say "candidates: realigning the $(grep -c '^>' "$OUT.cand.nets.fa") reads of the candidate families to $OUT.aug.fa"
  if ! "${RUSTLE_MINIMAP2:-minimap2}" -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t "$THREADS" "$OUT.aug.fa" "$OUT.cand.nets.fa" 2> "$OUT.aug.mm2.log" \
       | samtools sort -@ 2 -o "$OUT.aug.bam" -; then
    say "candidates: the patch realignment failed, see $OUT.aug.mm2.log"; rm -f "$OUT.aug.bam"; exit 1
  fi
  samtools index "$OUT.aug.bam"
  say "candidates: $OUT.aug.bam: $(samtools view -F 2308 "$OUT.aug.bam" | awk -F'\t' 'FILENAME == ARGV[1] { if (FNR > 1 && $5 == 1) c[$2] = 1; next } { n++; k += ($3 in c) } END { print n + 0 " primary records, " k + 0 " of them on a candidate contig" }' "$OUT.cand.candidates.tsv" -)"
}
# assign (O2): per-read copy assignment on the families' copy table (PREFIX.fam.copies.tsv/.fa, the default O1 output),
# swept over PREFIX.regions.txt (fam_regions). With --candidates and flagged candidates (PREFIX.aug.families.txt from the
# candidates stage of this families run): assign_with_candidates; candidate products present without --candidates are not
# used (one line says so). --legacy-catalog: the legacy catalog PREFIX.cat.* over whole contigs, as the stage ran before
# 2026-10-02. The three helpers are defined inside the stage (they are its code alone: figures/samples.py's driver_stage_code
# then keeps them out of the other stages' code hashes).
# o2 (register 1903): the family report of this run. Reads PREFIX.gtf (or the families' GTF) and PREFIX.fam.copies.tsv only.
stage_o2() {
  # RUSTLE_O2_EXTEND=1: the multi-copy extension first (near-tied multimappers seed the copies, floor-1 rescue inside families)
  # then the provenance test (register 1905 / 1906): a locus keeps its transcripts only when the EM over the tied reads leaves
  # it >= RUSTLE_O2_MIN_EM (2) reads -> PREFIX.o2.gtf (+ PREFIX.o2.em.tsv). The stage's own tie tables are not needed for it.
  # Default since 2026-11-22 (user: "lets readd it by default"); RUSTLE_O2_EXTEND=0 turns it off.
  # RUSTLE_O2_EM_WEIGHT (default 1, user 2026-11-23, register 1921): alignment scores inside the EM; 0 = counts-only. A failure of the extension
  # is reported and leaves the family report on the base assembly: it never fails the stage (or `all`).
  local ext=${RUSTLE_O2_EXTEND:-1}
  if [ "$ext" = 1 ]; then
    rm -f "$OUT.o2.gtf" "$OUT.o2.em.tsv"
    local tab=$OUT.molecules.tsv; [ -s "$tab" ] || tab=$BAM.as_table.tsv
    if ! ( RUSTLE_MULTICOPY_EQ=${RUSTLE_MULTICOPY_EQ:-0} stage_multicopy ); then
      say "o2: the multi-copy extension FAILED (see $OUT.multicopy.log); the report is on $OUT.gtf"; ext=0
    elif ! "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/../bench/o2_families/em_filter.py" "$OUT.multicopy.gtf" "$BAM" "$tab" "$OUT.o2.tmp" --min-em "${RUSTLE_O2_MIN_EM:-2}" $([ "${RUSTLE_O2_DROP:-0}" = 1 ] || echo --label) $([ "${RUSTLE_O2_EM_WEIGHT:-1}" = 0 ] || echo --weighted "${RUSTLE_O2_EM_WEIGHT:-1}") --ties "$OUT.o2.ties.tsv.gz" --threads "$JOB_THREADS" > "$OUT.o2.em.log" 2>&1; then
      say "o2: the provenance test FAILED (see $OUT.o2.em.log); the report is on $OUT.gtf"; rm -f "$OUT.o2.tmp".*; ext=0
    else
      mv -f "$OUT.o2.tmp.gtf" "$OUT.o2.gtf"; mv -f "$OUT.o2.tmp.em.tsv" "$OUT.o2.em.tsv"
      say "o2: $(tail -1 "$OUT.o2.em.log" | sed 's/ -> .*//') -> $OUT.o2.gtf"
      # The reads it collects are kept in PREFIX.o2.sequences.cache (RUSTLE_O2_SEQUENCES_CACHE=0 = no cache file): a rerun of
      # the tool on the same inputs with other options takes seconds to minutes (user 2026-11-23: keep runs re-runnable).
      # RUSTLE_O2_SEQUENCES=1 (opt-in, register 1926): IsoCon's test on the tied reads of each group of copies the reference
      # cannot separate -> PREFIX.o2.sequences.tsv and the column tied_sequences of the report. Never fails the stage.
      rm -f "$OUT.o2.sequences.tsv"
      if [ "${RUSTLE_O2_SEQUENCES:-0}" = 1 ]; then
        if [ ! -s "${FASTA:-}" ]; then say "o2: RUSTLE_O2_SEQUENCES needs --fasta; skipped"
        elif ! "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/../bench/o2_families/recurrence_test.py" "$OUT.multicopy.gtf" "$BAM" "$tab" "$FASTA" "$OUT.o2.sequences.tsv" --threads "$JOB_THREADS" $([ "${RUSTLE_O2_SEQUENCES_CACHE:-1}" = 0 ] || echo --cache "$OUT.o2.sequences.cache") --ties "$OUT.o2.ties.tsv.gz" > "$OUT.o2.sequences.log" 2>&1; then
          say "o2: the sequence test FAILED (see $OUT.o2.sequences.log); the report has no tied_sequences column"; rm -f "$OUT.o2.sequences.tsv" "$OUT.o2.sequences.tsv.tmp"
        else say "o2: $(tail -1 "$OUT.o2.sequences.log" | sed 's/ -> .*//') -> $OUT.o2.sequences.tsv"; fi
      fi
    fi
  fi
  local gtf=${FAM_GTF:-$OUT.gtf}; [ -s "$gtf" ] || gtf=$OUT.gtf
  [ -s "$gtf" ] && [ -s "$OUT.fam.copies.tsv" ] || { echo "[rustle_pipeline] o2 needs $OUT.gtf and $OUT.fam.copies.tsv (run assemble and families first)" >&2; exit 2; }
  local eqs=("$OUT.eq.loci.tsv" "$OUT.multicopy.eq.loci.tsv")
  if [ "$ext" = 1 ] && [ -s "$OUT.o2.gtf" ]; then gtf=$OUT.o2.gtf; eqs=("$OUT.multicopy.eq.loci.tsv" "$OUT.eq.loci.tsv"); fi
  local sqa=(); if [ "$ext" = 1 ] && [ -s "$OUT.o2.sequences.tsv" ]; then sqa=(--sequences "$OUT.o2.sequences.tsv"); fi
  local eq=(); local e; for e in "${eqs[@]}"; do
    if [ -s "$e" ] && [ ! "$e" -ot "$OUT.fam.copies.tsv" ]; then eq=(--eq "$e"); break; fi; done
  "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/o2_family_report.py" "$gtf" "$OUT.fam.copies.tsv" "$OUT.o2" "${eq[@]}" "${sqa[@]}" > "$OUT.o2.log" 2>&1 \
    || { echo "[rustle_pipeline] o2 failed (see $OUT.o2.log)" >&2; exit 1; }
  say "$(tail -1 "$OUT.o2.log") -> $OUT.o2.families.tsv, $OUT.o2.copies.tsv${eq:+ (with tie groups from ${eq[1]})}"
}
# o2bench (register 1903, opt-in): planted reads at the named families of the annotation, every assembler, scored against the truth.
stage_o2bench() {
  local b=$RP_DRIVER_DIR/../bench/o2_families d=$OUT.o2bench.d sc rc=0 scs=${RUSTLE_O2BENCH_SCENARIOS:-E,D,S}
  [ -s "$FASTA" ] || { echo "[rustle_pipeline] o2bench needs --fasta" >&2; exit 2; }
  mkdir -p "$d"
  if [ ! -s "$d/plant/plan.tsv" ]; then
    if [ -n "${GFF:-}" ]; then
      "${RUSTLE_PYTHON:-python3}" "$b/plant.py" "$GFF" "$FASTA" "$d/plant" --scenarios "$scs" --reads "${RUSTLE_O2BENCH_READS:-60}" > "$d/plant.log" 2>&1
    else
      [ -s "$OUT.gtf" ] && [ -s "$OUT.fam.copies.tsv" ] || { echo "[rustle_pipeline] o2bench without --gff plants at this run's own families: it needs $OUT.gtf and $OUT.fam.copies.tsv" >&2; exit 2; }
      "${RUSTLE_PYTHON:-python3}" "$b/plant.py" --from-run "$OUT.gtf" "$OUT.fam.copies.tsv" "$FASTA" "$d/plant" --scenarios "$scs" --reads "${RUSTLE_O2BENCH_READS:-30}" > "$d/plant.log" 2>&1
    fi || { echo "[rustle_pipeline] o2bench: planting failed (see $d/plant.log)" >&2; exit 1; }
  fi
  for sc in ${scs//,/ }; do
    # the child assemblies are runs of their own: none of this run's stage switches reach them
    env -u RP_ALL -u RP_LAYERS_BG RUSTLE_BIN="$BIN" bash "$b/run_tools.sh" "$d/plant" "$sc" "$FASTA" "$d/$sc" "$JOB_THREADS" > "$d/$sc.log" 2>&1 || rc=1
    [ -s "$d/$sc/score.tsv" ] || rc=1
  done
  [ "$rc" = 0 ] || { echo "[rustle_pipeline] o2bench: a scenario failed (see $d/*.log)" >&2; exit 1; }
  { head -1 "$d/${scs%%,*}/score.tsv"; for sc in ${scs//,/ }; do grep -P '\tALL\t' "$d/$sc/score.tsv"; done; for sc in ${scs//,/ }; do grep -v '^scenario' "$d/$sc/score.tsv" | grep -vP '\tALL\t'; done; } > "$OUT.o2bench.tsv"
  say "o2bench: $(grep -cP '\tALL\t' "$OUT.o2bench.tsv") tool x scenario rows -> $OUT.o2bench.tsv (planted truth in $d/plant)"
  grep -P '\tALL\t' "$OUT.o2bench.tsv" | awk -F'\t' '{printf "[rustle_pipeline]   %s %-18s copies %s/%s (recall %s, precision %s) isoforms recall %s precision %s\n", $1, $2, $5, $4, $6, $9, $12, $15}' >&2
}
stage_members() {
  local gtf=${FAM_GTF:-$OUT.gtf}; [ -s "$gtf" ] || gtf=$OUT.gtf
  [ -s "$gtf" ] || { echo "[rustle_pipeline] members needs $OUT.gtf (run assemble first)" >&2; exit 2; }
  # the run's score table gives each read's best AS (saves a pass); the satellite list assign used keeps the scan off the arrays
  local tab=$OUT.molecules.tsv; [ -s "$tab" ] || tab=$BAM.as_table.tsv; local mtab=(); [ ! -s "$tab" ] || mtab=(--as-table "$tab")
  local bed=${RUSTLE_EXCLUDE_BED-$FASTA.exclude.bed}; local mex=(); [ ! -s "$bed" ] || mex=(--exclude "$bed")
  "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/member_discovery.py" "$gtf" "$BAM" "$OUT.members.tmp" --min-reads "${RUSTLE_MEMBERS_MIN_READS:-5}" --min-ratio "${RUSTLE_MEMBERS_MIN_RATIO:-0.3}" "${mtab[@]}" "${mex[@]}" --threads "$JOB_THREADS" > "$OUT.members.log" 2>&1 \
    || { echo "[rustle_pipeline] members failed (see $OUT.members.log); no other product is affected" >&2; rm -f "$OUT.members.tmp".*; exit 1; }
  mv -f "$OUT.members.tmp.members.tsv" "$OUT.members.tsv"; mv -f "$OUT.members.tmp.members.gtf" "$OUT.members.gtf"
  say "$(tail -1 "$OUT.members.log" | sed 's/ -> .*//') -> $OUT.members.tsv, $OUT.members.gtf"
}
# provenance (opt-in, never in `all`; advisor's request 2026-11-23): where each read went. Reconstructed after the run from
# the BAM and the assembly (the assembler keeps no per-read record): tools/read_provenance.py.
stage_provenance() {
  local gtf=$OUT.o2.gtf; [ -s "$gtf" ] || gtf=$OUT.gtf
  [ -s "$gtf" ] || { echo "[rustle_pipeline] provenance needs $OUT.gtf (run assemble first)" >&2; exit 2; }
  local tab=$OUT.molecules.tsv; [ -s "$tab" ] || tab=$BAM.as_table.tsv; local a=(); [ ! -s "$tab" ] || a+=(--as-table "$tab")
  [ ! -s "${FASTA:-}" ] || a+=(--fasta "$FASTA"); [ ! -s "$OUT.shortexons.reads.tsv" ] || a+=(--shortexons "$OUT.shortexons.reads.tsv")
  [ ! -s "$OUT.o2.sequences.tsv" ] || a+=(--sequences "$OUT.o2.sequences.tsv"); [ ! -s "$OUT.fam.copies.tsv" ] || a+=(--copies "$OUT.fam.copies.tsv")
  [ ! -s "$OUT.polish_reasons.tsv" ] || a+=(--polish-reasons "$OUT.polish_reasons.tsv")
  "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/read_provenance.py" "$gtf" "$BAM" "$OUT.provenance" "${a[@]}" --threads "$JOB_THREADS" > "$OUT.provenance.log" 2>&1 \
    || { echo "[rustle_pipeline] provenance failed (see $OUT.provenance.log); no other product is affected" >&2; exit 1; }
  say "$(tail -1 "$OUT.provenance.log" | sed 's/ -> .*//') -> $OUT.provenance.reads.tsv.gz, $OUT.provenance.transcripts.tsv (on $gtf)"
}
stage_assign() {
  unset RUSTLE_READTHROUGH_JUNCTIONS   # assemble-only: copy_assign refuses it on every assign call (round-23 audit)
  # fam_regions all|only|skip [LIST]: copy_assign's --regions for the families of PREFIX.fam.copies.regions (`{fid}\t{chrom}:
  # {lo}-{hi}`, the copies' hull +- 5 kb per family and contig, as mcl_families writes it): every row, the rows of the
  # families LIST names (one id per line), or the other rows, MERGED per contig into disjoint intervals (an interval that
  # overlaps or touches the previous one joins it). copy_assign refuses overlapping --regions, and families' hulls
  # overlap and nest wherever families interleave (human chr16: 92 overlapping neighbours among 111 regions), so the
  # second column cannot be passed as it is. Every family's own interval lies inside exactly one merged region, which is
  # what binds it; its reads are still gathered around each copy (copy_assign's copy windows), not over the region.
  # `utilities candidate-augment` merges the candidate families' rows by the same rule (PREFIX.aug.regions.txt).
  fam_regions() {
    local mode=$1 list=${2:-/dev/null}
    awk -F'\t' -v OFS='\t' -v mode="$mode" -v list="$list" '
      FILENAME == list { if ($1 != "") pick[$1] = 1; next }
      mode == "all" || (mode == "only") == ($1 in pick) {
        p = match($2, /:[0-9]+-[0-9]+$/)
        if (!p) { print "[rustle_pipeline] bad region \"" $2 "\" in " FILENAME > "/dev/stderr"; bad = 1; exit 2 }
        split(substr($2, p + 1), r, "-"); print substr($2, 1, p - 1), r[1], r[2] }
      END { if (bad) exit 2 }' "$list" "$OUT.fam.copies.regions" |
      LC_ALL=C sort -t "$(printf '\t')" -k1,1 -k2,2n -k3,3n |
      awk -F'\t' '$1 != c || $2 > e { if (c != "") print c ":" s "-" e; c = $1; s = $2; e = $3; next }
                  $3 > e { e = $3 }
                  END { if (c != "") print c ":" s "-" e }'
  }
  # cand_ready: true when this run uses the candidates stage's products (--candidates, PREFIX.cand.candidates.tsv present);
  # a table older than PREFIX.fam.copies.tsv was made for other families: an error, not a silent mismatch.
  cand_ready() {
    [ "$CANDIDATES" = 1 ] && [ -s "$OUT.cand.candidates.tsv" ] || return 1
    if [ "$OUT.cand.candidates.tsv" -ot "$OUT.fam.copies.tsv" ]; then
      echo "[rustle_pipeline] $OUT.cand.candidates.tsv is older than $OUT.fam.copies.tsv (families ran again): run the candidates stage again, or drop --candidates" >&2; exit 2
    fi
  }
  # assign_with_candidates: the two-run O2 split (spec §4/§7; v1, §10 names the single run over a merged BAM as the end
  # state). The families PREFIX.aug.families.txt lists are assigned on the augmented inputs (PREFIX.aug.bam = their nets
  # realigned to PREFIX.aug.fa, PREFIX.aug.copies.* = their copies + candidates, PREFIX.aug.regions.txt), every other family
  # on the original BAM and copy table; copy_assign's --only-families / --skip-families select them, so each family is
  # assigned in exactly one run (a run with no family left is not made). Every per-family table (PREFIX.assign_cand.<t>.tsv
  # and PREFIX.assign_rest.<t>.tsv: assignments, families, quant, family_join, famcn_readonly, and with --inspect psv_reads,
  # psv_cols, psv_copies, posterior) is concatenated, one header, into PREFIX.assign.<t>.tsv; the run certificates are not
  # family tables and stay per run (PREFIX.assign_cand.params.tsv, PREFIX.assign_rest.params.tsv); PREFIX.assign.log holds
  # both logs. Ruling R13: each candidate family is assigned as a family that includes its candidate copies, under the
  # unchanged AS-tied gate, so a read the realignment placed uniquely on a candidate is not AS-tied and has no row (its
  # placement is PREFIX.aug.bam's); rows arise only for molecules tied between copies. How the split differs from one run
  # over every family: copy_assign keeps a process-wide registry of the molecules tied outside the catalog (keyed by read
  # name) and demotes a family's verdicts with it, so a row's `status`, not only its `tie_outside_catalog`, can differ from
  # a single run; the cross-family steps see only their own run's families (RUSTLE_XFAM_RECONCILE's reconciliation, and
  # `--union-certificate` in a hand-made split); and in the candidate run, candidate families whose real-chromosome hulls
  # overlap register each other's tied molecules and can demote each other's verdicts.
  assign_with_candidates() {
    [ -s "$OUT.aug.bam.bai" ] && [ ! "$OUT.aug.bam.bai" -ot "$OUT.aug.families.txt" ] \
      || { echo "[rustle_pipeline] $OUT.aug.bam is missing or older than $OUT.aug.families.txt: run the candidates stage again, or drop --candidates" >&2; exit 2; }
    local n_cand n_rest t f
    n_cand=$(grep -c . "$OUT.aug.families.txt")
    n_rest=$(awk -F'\t' 'FILENAME == ARGV[1] { c[$1] = 1; next } FNR > 1 && !($1 in c) { print $1 }' "$OUT.aug.families.txt" "$OUT.fam.copies.tsv" | sort -u | wc -l)
    say "assign: $n_cand candidate families on $OUT.aug.*, the other $n_rest on $OUT.fam.copies.tsv"
    rm -f "$OUT".assign.*.tsv "$OUT".assign_cand.*.tsv "$OUT".assign_rest.*.tsv
    "$BIN/copy_assign" "${SOLE[@]}" --bam "$OUT.aug.bam" --fasta "$OUT.aug.fa" --regions "$OUT.aug.regions.txt" --families "$OUT.aug.copies.tsv" \
      --copies-fa "$OUT.aug.copies.fa" --only-families "$OUT.aug.families.txt" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_cand" > "$OUT.assign_cand.log" 2>&1 \
      || { local rc=$?; say "assign: copy_assign failed on the candidate families (exit $rc), see $OUT.assign_cand.log; the likely cause: a real copy whose reads all realigned to its candidate, so it has no read in $OUT.aug.bam (copy_assign requires one on every real copy)"; exit "$rc"; }
    local runs=("$OUT.assign_cand")
    if [ "$n_rest" -gt 0 ]; then
      fam_regions skip "$OUT.aug.families.txt" > "$OUT.regions.txt"
      "$BIN/copy_assign" "${SOLE[@]}" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
        --copies-fa "$OUT.fam.copies.fa" --skip-families "$OUT.aug.families.txt" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_rest" > "$OUT.assign_rest.log" 2>&1 \
        || { local rc=$?; say "assign: copy_assign failed on the other families (exit $rc), see $OUT.assign_rest.log"; exit "$rc"; }
      runs+=("$OUT.assign_rest")
    fi
    for t in $(for r in "${runs[@]}"; do for f in "$r".*.tsv; do [ -e "$f" ] || continue; f=${f#"$r."}; echo "${f%.tsv}"; done; done | sort -u); do
      [ "$t" = params ] && continue
      local first=""
      : > "$OUT.assign.$t.tsv"
      for r in "${runs[@]}"; do
        [ -e "$r.$t.tsv" ] || continue
        if [ -z "$first" ]; then
          first=$r; cat "$r.$t.tsv" >> "$OUT.assign.$t.tsv"
        elif [ "$(head -1 "$r.$t.tsv")" = "$(head -1 "$first.$t.tsv")" ]; then
          awk 'NR > 1' "$r.$t.tsv" >> "$OUT.assign.$t.tsv"
        else
          echo "[rustle_pipeline] $r.$t.tsv and $first.$t.tsv have different headers: cannot concatenate them into $OUT.assign.$t.tsv" >&2
          rm -f "$OUT".assign.*.tsv; exit 2
        fi
      done
    done
    for r in "${runs[@]}"; do echo "=== [rustle_pipeline] $r.log"; cat "$r.log"; done > "$OUT.assign.log"
    # assigned_copy is the family's copy INDEX: family_join names the copy (copy_tid)
    local on_cand
    on_cand=$(awk -F'\t' 'FNR == 1 { next } FILENAME == ARGV[1] { tid[$1 "\t" $2] = $3; next } $4 == "assigned" && tid[$2 "\t" $3] ~ /^cand_/ { n++ } END { print n + 0 }' \
      "$OUT.assign.family_join.tsv" "$OUT.assign.assignments.tsv")
    say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family; $on_cand on a candidate copy; a read the realignment placed uniquely on a candidate is not AS-tied and has no row)"
  }
  if [ "$LEGACY_CATALOG" = 1 ]; then
    say "assign: per-read copy assignment on $OUT.cat"
    samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":0-"$3}' > "$OUT.regions.txt"
    "$BIN/copy_assign" "${SOLE[@]}" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
      --families "$OUT.cat.copies.tsv" --copies-fa "$OUT.cat.copies.fa" "${ASSIGN_FAST[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1
    say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
    return 0
  fi
  [ -s "$OUT.fam.copies.tsv" ] || { echo "assign needs $OUT.fam.copies.tsv (run the families stage with an mcl_families that writes the copy table), or pass --legacy-catalog" >&2; exit 2; }
  if [ "$(awk 'NR > 1 && NF' "$OUT.fam.copies.tsv" | wc -l)" = 0 ]; then
    rm -f "$OUT".assign.*.tsv "$OUT".assign.*.txt "$OUT".assign_cand.* "$OUT".assign_rest.*   # an earlier run's products
    say "assign: $OUT.fam.copies.tsv lists no copy (no family): nothing to assign"; return 0
  fi
  [ -s "$OUT.fam.copies.regions" ] || { echo "assign needs $OUT.fam.copies.regions (the families stage writes it beside the copy table)" >&2; exit 2; }
  # flagged candidates but no augmented families: the candidates stage did not finish; assigning without them would be silent
  if cand_ready && [ ! -s "$OUT.aug.families.txt" ] && awk -F'\t' 'NR > 1 && $5 == 1 { f = 1 } END { exit !f }' "$OUT.cand.candidates.tsv"; then
    echo "[rustle_pipeline] $OUT.cand.candidates.tsv has flagged candidates but $OUT.aug.families.txt is missing (the candidates stage did not finish): run it again, or drop --candidates" >&2; exit 2
  fi
  if cand_ready && [ -s "$OUT.aug.families.txt" ]; then
    # the candidates path applies none of the skip lists below: refused rather than silently ignored (round-10 audit)
    if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}${RUSTLE_EXCLUDE_SOFTMASK:-}${RUSTLE_EXCLUDE_FASTA:-}" ]; then
      echo "[rustle_pipeline] --candidates does not apply --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED / RUSTLE_EXCLUDE_SOFTMASK / RUSTLE_EXCLUDE_FASTA: drop one of them" >&2; exit 2
    fi
    assign_with_candidates; return 0
  fi
  if [ "$CANDIDATES" != 1 ] && [ -s "$OUT.aug.families.txt" ] && [ ! "$OUT.aug.families.txt" -ot "$OUT.fam.copies.tsv" ]; then
    say "assign: the candidates stage's products ($OUT.aug.*) are present but unused (the stage is opt-in, ruling R14): pass --candidates to assign the candidate families on them"
  fi
  say "assign: per-read copy assignment on $OUT.fam.copies.tsv"
  rm -f "$OUT".assign.*.tsv "$OUT".assign_cand.* "$OUT".assign_rest.*   # an earlier split run's tables
  local mem=()
  # --tied-prefilter (byte-identical): each region is read light first and only the AS-tied molecules' records are loaded in full
  mem+=(--tied-prefilter)
  # --max-window-records N (--low-mem: 500,000): a family whose merged copy windows (copy +- 50 kb clipped to its region, what copy_assign reads) hold more
  # than N records copy_assign would keep (rDNA arrays, acrocentric short arms) is not assigned, and its regions are not swept: listed
  # with the count reached in PREFIX.assign.depth_skipped.tsv (`utilities window-depth`). Only reads such a family shares with others
  # can change: their placement on its copies is no longer seen.
  rm -f "$OUT.assign.depth_skipped.tsv" "$OUT.assign.depth_skipped.txt"
  if [ "$MAX_WINDOW_RECORDS" -gt 0 ]; then
    # clipped to the regions copy_assign sweeps (every family's: the skip lists only shrink them, so still an upper bound)
    fam_regions all > "$OUT.assign.depth_regions.txt"
    "$BIN/utilities" window-depth --bam "$BAM" --copies "$OUT.fam.copies.tsv" --max "$MAX_WINDOW_RECORDS" \
      --regions "$OUT.assign.depth_regions.txt" --out "$OUT.assign.depth_skipped.tsv" 2>> "$OUT.assign.depth.log"
    if [ -s "$OUT.assign.depth_skipped.tsv" ]; then
      cut -f1 "$OUT.assign.depth_skipped.tsv" > "$OUT.assign.depth_skipped.txt"   # --skip-families is passed once below
      say "assign: $(wc -l < "$OUT.assign.depth_skipped.tsv") famil(ies) with a copy window over $MAX_WINDOW_RECORDS records are not assigned ($OUT.assign.depth_skipped.tsv)"
    fi
  fi
  # RUSTLE_ANNOTATED_ONLY=GFF (opt-in, row 1551): families with no copy on a `gene` feature of GFF join the skip list
  # (A119b: 212 of 1,565, among them every chr13 hotspot family that sets the assign wall). Output-changing: not a default.
  rm -f "$OUT.assign.unannotated.txt"
  if [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ]; then
    "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/annotated_families.py" "$OUT.fam.copies.tsv" "$RUSTLE_ANNOTATED_ONLY" > "$OUT.assign.unannotated.txt" 2>> "$OUT.assign.depth.log" || { echo "[rustle_pipeline] assign: the skip list could not be built (see $OUT.assign.depth.log)" >&2; exit 2; }
    say "assign: annotated-only: $(wc -l < "$OUT.assign.unannotated.txt") famil(ies) with no copy on a gene of $RUSTLE_ANNOTATED_ONLY are not assigned ($OUT.assign.unannotated.txt)"
  fi
  # RUSTLE_EXCLUDE_BED=BED (opt-in, row 1559): families whose copies mostly lie in a confounding element (T2T censat
  # satellites: the chr13 HSat1A array is 618 of A119b's 1,565 families and the assign wall) join the skip list.
  # Use a SATELLITE-ONLY list (user decision 2026-11-15, rows 1749 / 1750): the whole censat v2.1 BED also holds `ct`
  # transition rows (NBPF, POTE, GOLGA) and generic `censat` rows (the TSPY / RBMY arrays), and the families it names held
  # 682 of A119b's 1,765 true pairs. docs/DATA.md has the list (hsat / hor / dhor / mon / bsat / gsat / rDNA rows) and its command.
  # The genome's own list as a sidecar, FASTA.exclude.bed, is used when RUSTLE_EXCLUDE_BED is unset (user 2026-11-15, as the
  # SEDEF sidecar of `families`); RUSTLE_EXCLUDE_BED= (empty) turns it off. CHM13's satellite-only list is linked that way.
  if [ -z "${RUSTLE_EXCLUDE_BED+x}" ] && [ -s "$FASTA.exclude.bed" ]; then RUSTLE_EXCLUDE_BED=$FASTA.exclude.bed; fi
  rm -f "$OUT.assign.excluded.txt"
  if [ -n "${RUSTLE_EXCLUDE_BED:-}" ]; then
    if awk -F'\t' '$4 ~ /^ct_/ {f=1; exit} END{exit !f}' "$RUSTLE_EXCLUDE_BED" 2>/dev/null; then
      say "assign: WARNING $RUSTLE_EXCLUDE_BED has centromeric-transition (ct_) rows: gene families there (NBPF, POTE, GOLGA) will not be assigned; use the satellite-only list (docs/DATA.md)"
    fi
    # judged on the first run's copies when famreps added representatives: a family must not fall under the rule because of
    # them (row 1662: MCL215, 8 -> 26 copies, lost its 129 assigned rows that way)
    local excl_tbl="$OUT.fam.copies.tsv"; [ -s "$OUT.fam.base.copies.tsv" ] && excl_tbl="$OUT.fam.base.copies.tsv"
    "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/annotated_families.py" "$excl_tbl" "$RUSTLE_EXCLUDE_BED" --majority-in > "$OUT.assign.excluded.txt" 2>> "$OUT.assign.depth.log" || { echo "[rustle_pipeline] assign: the skip list could not be built (see $OUT.assign.depth.log)" >&2; exit 2; }
    say "assign: exclude-bed: $(wc -l < "$OUT.assign.excluded.txt") famil(ies) mostly inside $RUSTLE_EXCLUDE_BED are not assigned ($OUT.assign.excluded.txt)"
  fi
  # RUSTLE_EXCLUDE_SOFTMASK=FRAC (opt-in, row 1690): the same hatch from the FASTA's soft mask instead of a BED: families
  # with most copies at >= FRAC lowercase exon bases. Judged on the first run's copies, as the BED rule is.
  rm -f "$OUT.assign.softmask.txt"
  if [ -n "${RUSTLE_EXCLUDE_SOFTMASK:-}" ]; then
    local sm_tbl="$OUT.fam.copies.tsv"; [ -s "$OUT.fam.base.copies.tsv" ] && sm_tbl="$OUT.fam.base.copies.tsv"
    "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/softmask_families.py" "$sm_tbl" "$FASTA" "$RUSTLE_EXCLUDE_SOFTMASK" > "$OUT.assign.softmask.txt" 2>> "$OUT.assign.depth.log" || { echo "[rustle_pipeline] assign: the soft-mask skip list could not be built (see $OUT.assign.depth.log)" >&2; exit 2; }
    say "assign: exclude-softmask: $(wc -l < "$OUT.assign.softmask.txt") famil(ies) with most copies >= $RUSTLE_EXCLUDE_SOFTMASK soft-masked are not assigned ($OUT.assign.softmask.txt)"
  fi
  # RUSTLE_EXCLUDE_FASTA=REF (opt-in, row 1691): the same hatch by sequence: families with most copies >= 50% covered by
  # alignments to REF (an rDNA unit). Judged on the first run's copy sequences, as the other two rules are.
  rm -f "$OUT.assign.seqexcl.txt"
  if [ -n "${RUSTLE_EXCLUDE_FASTA:-}" ]; then
    [ -s "$RUSTLE_EXCLUDE_FASTA" ] || { echo "[rustle_pipeline] assign: RUSTLE_EXCLUDE_FASTA=$RUSTLE_EXCLUDE_FASTA is missing or empty" >&2; exit 2; }
    local sq_fa="$OUT.fam.copies.fa"; [ -s "$OUT.fam.base.copies.fa" ] && sq_fa="$OUT.fam.base.copies.fa"
    RUSTLE_SEQFAM_THREADS="$THREADS" "${RUSTLE_PYTHON:-python3}" "$RP_DRIVER_DIR/sequence_families.py" "$sq_fa" "$RUSTLE_EXCLUDE_FASTA" > "$OUT.assign.seqexcl.txt" 2>> "$OUT.assign.depth.log" || { echo "[rustle_pipeline] assign: the sequence skip list could not be built (see $OUT.assign.depth.log)" >&2; exit 2; }
    say "assign: exclude-fasta: $(wc -l < "$OUT.assign.seqexcl.txt") famil(ies) with most copies covered by $RUSTLE_EXCLUDE_FASTA are not assigned ($OUT.assign.seqexcl.txt)"
  fi
  # the regions are rebuilt without the skipped families: a region left with no family would be read whole
  { for f in "$OUT.assign.depth_skipped.txt" "$OUT.assign.unannotated.txt" "$OUT.assign.excluded.txt" "$OUT.assign.softmask.txt" "$OUT.assign.seqexcl.txt"; do [ -s "$f" ] && cat "$f"; done; true; } | sort -u > "$OUT.assign.skipped.txt"
  if [ -s "$OUT.assign.skipped.txt" ]; then
    mem+=(--skip-families "$OUT.assign.skipped.txt")
    fam_regions skip "$OUT.assign.skipped.txt" > "$OUT.regions.txt"
  else
    fam_regions all > "$OUT.regions.txt"
  fi
  if [ ! -s "$OUT.regions.txt" ]; then
    say "assign: every family is on a skip list ($OUT.assign.skipped.txt): nothing to assign"; return 0
  fi
  "$BIN/copy_assign" "${SOLE[@]}" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" --families "$OUT.fam.copies.tsv" \
    --copies-fa "$OUT.fam.copies.fa" "${ASSIGN_FAST[@]}" "${mem[@]}" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign" > "$OUT.assign.log" 2>&1 \
    || { local rc=$?; say "assign: copy_assign failed (exit $rc), see $OUT.assign.log"; exit "$rc"; }
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign.assignments.tsv" | wc -l) (one row per read x family)"
}
# assign on the FAMILIES copy table (the merged stage's phase 3): the same byte-identical levers the genome-wide
# sweep pre-registers recommended — `--skip-poa-diagnostic` (the POA homology pass is purely diagnostic; measured
# 6.8x on the heaviest family) and `--region-threads` (cross-contig rayon pool; byte-identical by construction) —
# are ON here by default instead of opt-in. Products: PREFIX.assign_fam.* (the catalog assign keeps PREFIX.assign.*).
stage_assign_fam() {
  # this phase reads every family: the skip lists of the `assign` stage are not applied, so they are refused (round-10 audit)
  if [ "$MAX_WINDOW_RECORDS" -gt 0 ] || [ -n "${RUSTLE_ANNOTATED_ONLY:-}" ] || [ -n "${RUSTLE_EXCLUDE_BED:-}${RUSTLE_EXCLUDE_SOFTMASK:-}${RUSTLE_EXCLUDE_FASTA:-}" ]; then
    echo "[rustle_pipeline] merged applies no --max-window-records / --low-mem / RUSTLE_ANNOTATED_ONLY / RUSTLE_EXCLUDE_BED / RUSTLE_EXCLUDE_SOFTMASK / RUSTLE_EXCLUDE_FASTA: use \`all\` for them" >&2; exit 2
  fi
  say "assign: per-read copy assignment on $OUT.fam.copies (families table; --skip-poa-diagnostic --region-threads $THREADS)"
  samtools view -H "$BAM" | awk '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2":0-"$3}' > "$OUT.regions.txt"
  # RUSTLE_READTHROUGH_JUNCTIONS is assemble-only (copy_assign refuses it here): dropped for this call alone, so the
  # merged.env record written after this phase still names it and a rerun resumes (round-24 audit)
  env -u RUSTLE_READTHROUGH_JUNCTIONS "$BIN/copy_assign" "${SOLE[@]}" --bam "$BAM" --fasta "$FASTA" --regions "$OUT.regions.txt" \
    --families "$OUT.fam.copies.tsv" --copies-fa "$OUT.fam.copies.fa" \
    --skip-poa-diagnostic --region-threads "$THREADS" "${INSPECT_ASSIGN[@]}" --out "$OUT.assign_fam" > "$OUT.assign_fam.log" 2>&1
  say "assign: $(awk -F'\t' 'NR>1 && $4=="assigned"' "$OUT.assign_fam.assignments.tsv" | wc -l) assigned rows of $(awk 'NR>1' "$OUT.assign_fam.assignments.tsv" | wc -l) (one row per read x family)"
}
# The genome-wide entry point: assemble -> families -> assign(fam), one command. Resume: PREFIX.merged.env records
# the exact inputs/settings of a COMPLETED merged run (written only when all three phases finished); a phase is
# skipped only when that file matches AND the phase's product exists, so a stale product from another BAM/genome is
# never silently reused. A crashed run redoes the phases (the families all-vs-all replays from PREFIX.cache when
# the loci FASTA is unchanged — the expensive step is never repeated).
# seed_sec= keeps its two values for the default and primary pools (a finished run of either is still recognised); the all pool and
# a non-default width of the good pool add a line, so a finished run of one pool never stands in for another
merged_env() {
  printf 'bam=%s\nfasta=%s\nthreads=%s\nseed_sec=%s\nbridge=%s\n' "$(readlink -f "$BAM")" "$(readlink -f "$FASTA")" "$THREADS" \
    "$([ "$SEED_POOL" = primary ] && echo 0 || echo 1)" "$BRIDGE_MODE"
  if [ "$SEED_POOL" = all ]; then printf 'seed_pool=all\n'; fi
  if [ "$SEED_POOL" = good ] && [ "$SEED_RATIO" != 0.98 ]; then printf 'seed_ratio=%s\n' "$SEED_RATIO"; fi
  if [ -n "$CONTIG" ]; then printf 'contig=%s\n' "$CONTIG"; fi
  # the SEDEF sidecar the families stage picks up when RUSTLE_SEDEF_BED is unset (round-4 audit)
  if [ -z "${RUSTLE_SEDEF_BED+x}" ] && [ -s "$FASTA.sedef.bed" ]; then stat -L -c "sedef_sidecar %s %Y" "$FASTA.sedef.bed"; fi
  if [ -z "${RUSTLE_EXCLUDE_BED+x}" ] && [ -s "$FASTA.exclude.bed" ]; then stat -L -c "exclude_sidecar %s %Y" "$FASTA.exclude.bed"; fi
  # every RUSTLE_* knob changes some stage's products: a resume must not skip a stage run with other values (audit 2026-11-07)
  # minus the speed knobs the driver derives from the allocation (threads, read budget): products are byte-identical across
  # them, and a resume in a differently sized job must not redo every stage (audit 2026-11-19)
  env | grep '^RUSTLE_' | grep -vE '^RUSTLE_(BAM_THREADS|REGION_READ_BUDGET|FAMILY_THREADS|STAR_THREADS|REGION_THREADS)=' | LC_ALL=C sort || true
  # and the binaries: a rebuild must rerun the stages it changes, not replay the old products (audit 2026-11-07, second pass)
  local b; for b in copy_assign mcl_families as_table utilities missing_copy_flag candidate_copies gw_family_catalog; do
    [ -n "${BIN:-}" ] && [ -e "$BIN/$b" ] && stat -L -c "bin $b %s %Y" "$BIN/$b"
  done; true
}
stage_merged() {
  # computed once, before stage_families sets RUSTLE_SEDEF_BED from the sidecar: the record written at the end used to lack
  # the sedef_sidecar line the next run's check has, so a genome with FASTA.sedef.bed never resumed
  local env_ok=0 env_now rec_t=0; env_now=$(merged_env)
  if [ -f "$OUT.merged.env" ] && [ "$(cat "$OUT.merged.env")" = "$env_now" ]; then env_ok=1; rec_t=$(stat -c %Y "$OUT.merged.env"); fi
  # a product rewritten after the record (a separate `assemble` / `families` on this prefix) is not the merged run's (audit 2026-11-19)
  recorded() { [ -s "$1" ] && [ "$(stat -c %Y "$1")" -le "$rec_t" ]; }
  # removed before any phase runs (a run that dies midway must not leave the record vouching for its products), and once one
  # phase runs every later phase runs too (its input changed), as in stage_all
  rm -f "$OUT.merged.env"
  if [ "$env_ok" = 1 ] && recorded "$OUT.gtf"; then
    say "merged: phase 1 assemble SKIPPED ($OUT.gtf present, env matches)"
  else
    stage_assemble; env_ok=0
  fi
  if [ "$env_ok" = 1 ] && recorded "$OUT.fam.copies.tsv"; then
    say "merged: phase 2 families SKIPPED ($OUT.fam.copies.tsv present, env matches)"
  else
    stage_families; env_ok=0
  fi
  if [ "$env_ok" = 1 ] && recorded "$OUT.assign_fam.assignments.tsv"; then
    say "merged: phase 3 assign SKIPPED ($OUT.assign_fam.assignments.tsv present, env matches)"
  else
    stage_assign_fam
  fi
  printf '%s\n' "$env_now" > "$OUT.merged.env"
  say "merged: done (env recorded in $OUT.merged.env; phases skip on the next identical run)"
}
# flag's first step alone: the loci's reads, mixture test and consensus (PREFIX.flag_scan.*). It reads only the BAM and the
# loci, an assemble product, so `all` runs it in the background from the end of assemble (round-39 audit: A119b families +
# assign + scan 890 s overlapped vs ~1,137 s serial, products identical).
stage_flag_scan() {
  [ -n "$GFF" ] || fam_gtf_guard flag
  local LOCI=${GFF:-$FAM_GTF}
  say "flag: scan $BAM on $LOCI"
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" ${GFF:+--gff "$GFF"} --index x --threads "$JOB_THREADS" \
    --out "$OUT.flag_scan" --scan-only > "$OUT.flag_scan.log" 2>&1
}
stage_flag() {
  [ -n "$INDEX" ] || { echo "flag needs --index (splice .mmi of the primary genome)" >&2; exit 2; }
  # the loci to scan: the annotation with --gff, else the de novo loci `families` reads (bridges are relations, not loci)
  [ -n "$GFF" ] || fam_gtf_guard flag
  local LOCI=${GFF:-$FAM_GTF} cand=()
  # with --candidates, the two O3 sources corroborate: the verdict table's `o3_candidate` column names a flagged candidate on
  # the row's locus (a candidates table older than PREFIX.fam.copies.tsv was made for other families: an error)
  if [ "$CANDIDATES" = 1 ] && [ -s "$OUT.cand.candidates.tsv" ]; then
    [ ! "$OUT.cand.candidates.tsv" -ot "$OUT.fam.copies.tsv" ] \
      || { echo "[rustle_pipeline] $OUT.cand.candidates.tsv is older than $OUT.fam.copies.tsv (families ran again): run the candidates stage again, or drop --candidates" >&2; exit 2; }
    cand=(--candidates "$OUT.cand.candidates.tsv")
  fi
  if [ "${RP_FLAG_SCAN_DONE:-0}" = 1 ]; then say "flag: scan ran in the background ($OUT.flag_scan.scan.tsv)"; else stage_flag_scan; fi
  say "flag: align + verdict"
  # a FASTA as --index (the --low-mem profile) is indexed on the fly in INDEX_GB-Gbp parts; a prebuilt .mmi is used as it is
  local igb=()
  case "$INDEX" in *.mmi) ;; *) igb=(--index-gb "$INDEX_GB");; esac
  "$BIN/missing_copy_flag" --bam "$BAM" --fasta "$FASTA" --loci "$LOCI" --index "$INDEX" --threads "$MM_THREADS" "${igb[@]}" \
    "${CONFIRM[@]}" "${FOREIGN[@]}" "${cand[@]}" --out "$OUT.flag" --from-scan "$OUT.flag_scan" > "$OUT.flag.log" 2>&1
  say "flag: $(grep -o 'verdicts: .*' "$OUT.flag.log")"
  if [ ${#cand[@]} -gt 0 ]; then say "flag: $(grep -o '[0-9]* flagged candidates with a nearest locus; [0-9]* rows name one' "$OUT.flag.log") (o3_candidate column)"; fi
}
# index: the splice index flag aligns its consensus sequences against (shared next to the FASTA, below). `-I INDEX_GB` Gbp per part: the default
# 8 gives one part for a mammalian genome (the figures' index, ~13 GB to hold); --index-gb 1 / --low-mem writes ~1-Gbp parts
# that flag loads one at a time (minimap2 --split-prefix merges their hits: near-identical, not byte-identical).
# The index depends on the FASTA, the part size and minimap2 alone, so it is SHARED by every run on that genome when it can be
# (round-2 audit: each run directory built its own ~13 GB copy): in RUSTLE_INDEX_DIR when set, else next to the FASTA unless
# the FASTA lies under $HOME (large files stay off the home quota), as <FASTA name>.splice.I<GB>G.mmi; otherwise, or when that
# directory takes no new file, PREFIX.splice.mmi. Its .key (FASTA path, size, mtime, part size, minimap2 version) must match
# to be reused; it is built under a temporary name and renamed, so a concurrent run never reads half of it.
index_key() { printf 'fasta=%s\nfasta_stat=%s\nindex_gb=%s\nminimap2=%s\n' "$(readlink -f "$FASTA")" "$(stat -Lc '%s %Y' "$FASTA")" "$INDEX_GB" "$("${RUSTLE_MINIMAP2:-minimap2}" --version 2>/dev/null)"; }
index_path() {
  local f d p probe; f=$(readlink -f "$FASTA"); d=${RUSTLE_INDEX_DIR:-$(dirname "$f")}
  if [ -z "${RUSTLE_INDEX_DIR:-}" ] && [[ "$f" == "$(readlink -f "$HOME")"/* ]]; then echo "$OUT.splice.mmi"; return; fi
  p="$d/$(basename "$f").splice.I${INDEX_GB}G.mmi"; probe="$p.probe.$$"
  if { [ -s "$p" ] && [ "$(cat "$p.key" 2>/dev/null)" = "$(index_key)" ]; } || { : > "$probe"; } 2>/dev/null; then rm -f "$probe"; echo "$p"
  else echo "$OUT.splice.mmi"; fi
}
stage_index() {
  local idx=${INDEX:-$(index_path)} key; key=$(index_key)
  if [ -s "$idx" ] && [ "$(cat "$idx.key" 2>/dev/null)" = "$key" ]; then say "index: $idx present"; return 0; fi
  rm -f "$idx.key"
  say "index: minimap2 -x splice -I ${INDEX_GB}G -d $idx"
  # a failed or killed build must not leave a ~13 GB partial index beside a shared genome (audit 2026-11-19)
  # added to the EXIT trap, not replacing it: the owner's trap removes the driver snapshot (line ~289)
  local prev_trap; prev_trap=$(trap -p EXIT | sed "s/^trap -- '//; s/' EXIT\$//")
  trap "rm -f '$idx.tmp.$$'; $prev_trap" EXIT
  "${RUSTLE_MINIMAP2:-minimap2}" -x splice -I "${INDEX_GB}G" -t "$JOB_THREADS" -d "$idx.tmp.$$" "$FASTA" > "$OUT.index.log" 2>&1 \
    || { rm -f "$idx.tmp.$$"; echo "[rustle_pipeline] index: minimap2 failed (see $OUT.index.log)" >&2; exit 1; }
  mv -f "$idx.tmp.$$" "$idx"; if [ -n "$prev_trap" ]; then trap "$prev_trap" EXIT; else trap - EXIT; fi
  printf '%s\n' "$key" > "$idx.key"
}
# stage_all_env STAGE: what a finished `all` stage was run with; a stage is skipped only when it matches (and its product is there)
stage_all_env() {
  merged_env
  printf 'candidates=%s\nlegacy=%s\ngff=%s\nindex=%s\nindex_gb=%s\nregion_threads=%s\npoa=%s\n' "$CANDIDATES" "$LEGACY_CATALOG" \
    "$GFF" "${INDEX:-$OUT.splice.mmi}" "$INDEX_GB" "${RUSTLE_REGION_THREADS:-8}" "${RUSTLE_POA_DIAGNOSTIC:-0}"
  printf 'max_window_records=%s\n' "$MAX_WINDOW_RECORDS"
  # settings only one stage reads go only into that stage's record (a flag-only option re-runs flag alone; round-23 audit),
  # and the driver itself everywhere (a default flipped here reaches no binary stat; round-10 audit)
  case "$1" in
    flag) printf 'confirm=%s\nforeign=%s\n' "${CONFIRM[*]-}" "${FOREIGN[*]-}";;
    candidates) printf 'delta=%s\ncand_max=%s\n' "$DELTA" "$CAND_MAX";;
    catalog|assign) printf 'inspect=%s\n' "$INSPECT";;
  esac
  [ -f "$0" ] && stat -L -c "driver %s %Y" "$0"; true
}
# stage_all: every stage as its own driver call under /usr/bin/time (PREFIX.<stage>.time: wall seconds, peak RSS of the stage's
# process tree, exit). Resume: a stage whose product exists and is not newer than its time record (a direct stage call rewrote
# it), whose time record says exit 0 and whose OWN settings record PREFIX.<stage>.env matches is skipped — until one stage runs,
# after which every later stage runs too. A stage's record is removed before it runs and written only when it succeeds, so a
# failed or killed stage is redone and the stages before it are not; `--from` writes no record for the stages it skips.
stage_all() {
  local plan=(assemble:"$OUT.gtf" families:"$OUT.fam.copies.tsv")
  if [ "$LEGACY_CATALOG" = 1 ]; then plan+=(catalog:"$OUT.cat.copies.tsv")
  elif [ "$CANDIDATES" = 1 ]; then plan+=(candidates:"$OUT.cand.candidates.tsv"); fi
  # register 1903 (user, 2026-11-22): the O2 stage of `all` is the family report; assign only on request or where its inputs are asked for
  # `all` never runs the legacy `assign` (no switch brings it back): after --candidates or --legacy-catalog it is called by name
  plan+=(o2:"$OUT.o2.families.tsv")
  local env_now ran=0 st prod args=("${ORIG_ARGS[@]}")
  # audit 2026-11-19: the child stages get this run's pool explicitly; without it a child `assemble` took the assemble-alone
  # default (primary, register 1812) while PREFIX.assemble.env recorded good
  if [ -n "${SEED_POOL:-}" ]; then args+=(--seed-pool "$SEED_POOL"); fi
  # GNU time: the system's, else one on PATH or pixi's (`pixi global install time`; compute nodes have no /usr/bin/time)
  local gtime=/usr/bin/time; [ -x "$gtime" ] || gtime=$(PATH="$PATH:${HOME:-}/.pixi/bin" type -P time || true)
  export RP_ALL=1   # the stage calls below are part of `all`: assemble leaves the guided assembly to this loop
  # families, assign, index and flag all run minimap2: a missing one fails at start, not after the assembly (round-13 audit)
  command -v "${RUSTLE_MINIMAP2:-minimap2}" > /dev/null || { echo "[rustle_pipeline] no ${RUSTLE_MINIMAP2:-minimap2} on PATH (families, assign and flag run it); set RUSTLE_MINIMAP2 or PATH" >&2; exit 2; }
  # round 52: so do assign's skip lists need python3 (tools/annotated_families.py); without one, assemble and families ran first
  # row 1662: so do families' extra representatives (inline, as stage_all is also run on its own by the dispatch tests)
  [ "${RUSTLE_FAMREPS:-1}" = 0 ] || command -v "${RUSTLE_PYTHON:-python3}" > /dev/null \
    || { echo "[rustle_pipeline] no ${RUSTLE_PYTHON:-python3} on PATH (families adds extra locus representatives with it); set RUSTLE_PYTHON or PATH, or RUSTLE_FAMREPS=0" >&2; exit 2; }
  if [ "$CANDIDATES" = 0 ] && [ -z "${RUSTLE_EXCLUDE_BED+x}" ] && [ -s "$FASTA.exclude.bed" ]; then export RUSTLE_EXCLUDE_BED=$FASTA.exclude.bed; fi   # the sidecar, as assign resolves it
  local v; for v in RUSTLE_EXCLUDE_BED RUSTLE_EXCLUDE_SOFTMASK RUSTLE_EXCLUDE_FASTA RUSTLE_ANNOTATED_ONLY; do
    if [ -n "${!v:-}" ] && ! command -v "${RUSTLE_PYTHON:-python3}" > /dev/null; then
      echo "[rustle_pipeline] no ${RUSTLE_PYTHON:-python3} on PATH ($v needs it for assign's skip list); set RUSTLE_PYTHON or PATH, or unset $v" >&2; exit 2
    fi
  done
  if [ "$NO_FLAG" = 0 ]; then
    # no --index: the default profile builds the shared FASTA.splice.I<GB>G.mmi (index stage); --low-mem aligns against the FASTA itself,
    # indexed on the fly in INDEX_GB-Gbp parts (no index file on disk)
    if [ -z "$INDEX" ] && [ "$LOW_MEM" = 1 ]; then INDEX=$FASTA; args+=(--index "$FASTA")
    elif [ -z "$INDEX" ]; then INDEX=$(index_path); plan+=(index:"$INDEX"); args+=(--index "$INDEX"); fi
    plan+=(flag:"$OUT.flag.missing_copy.tsv")
  fi
  local v; for v in "$FROM_STAGE" "$TO_STAGE"; do   # a misspelled --to ran every stage, a misspelled --from none
    [ -z "$v" ] || [[ " ${plan[*]%%:*} " == *" $v "* ]] || { echo "[rustle_pipeline] --from/--to '$v' is not a stage of this run (${plan[*]%%:*})" >&2; exit 2; }
  done
  local started=1; [ -n "$FROM_STAGE" ] && started=0
  # round-38 audit: the index depends on the FASTA alone, so a missing or stale one is built in the background from the start
  # (under assemble, which is BAM-bound) instead of after assign: ~130 s off a new genome's first run. Its step below waits.
  local idx_pid="" scan_pid="" flag_bg=0
  if [[ " ${plan[*]} " == *" index:"* ]] && [ -z "$FROM_STAGE$TO_STAGE" ] \
     && ! { [ -s "$INDEX" ] && [ "$(cat "$INDEX.key" 2>/dev/null)" = "$(index_key)" ]; }; then
    (
      if [ -n "$gtime" ]; then "$gtime" -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$OUT.index.time.tmp" bash "$0" index "${args[@]}"
      else st0=$SECONDS; rc=0; bash "$0" index "${args[@]}" || rc=$?; printf 'wall_s\t%s\n' "$(( SECONDS - st0 ))" > "$OUT.index.time.tmp"; exit "$rc"; fi
    ) 2> "$OUT.index.bg.err" &
    idx_pid=$!
    say "all: index building in the background ($INDEX)"
  fi
  # round 46: the guided assembly is a side product of assemble, written once assemble is behind us (run or skipped) and
  # outside the plan, so writing it never makes a later stage rerun
  local guided_due=0; [ -n "${GUIDE:-}" ] && guided_due=1
  # round 60: a process and everything below it, children first found then killed (flag's align has a minimap2 child that
  # `pkill -P` on the background driver did not reach)
  kill_tree() { local c; for c in $(pgrep -P "$1" 2>/dev/null); do kill_tree "$c"; done; kill "$1" 2>/dev/null || true; }
  # round 66: the layered assembly finishing in the background (PREFIX.final.pid, not a child of this shell)
  [ "$LOW_MEM" = 1 ] || export RP_LAYERS_BG=1
  # the layered run of PREFIX.final.pid is alive only when its pid runs on this host with the recorded start time
  final_alive() {
    local p h t; read -r p h t < "$OUT.final.pid" 2>/dev/null || return 1
    [ -n "$p" ] && [ "$h" = "$(hostname)" ] && [ "$(ps -o lstart= -p "$p" 2>/dev/null | tr -s ' ' '_')" = "$t" ] && kill -0 "$p" 2>/dev/null
  }
  all_layers_wait() {
    [ -s "$OUT.final.pid" ] || return 0
    while final_alive; do sleep 0.2; done
    if [ ! -s "$OUT.final.rc" ]; then
      # it ended without writing its code (killed): the assembly has to run again
      rm -f "$OUT.final.pid" "$OUT.assemble.env"
      echo "[rustle_pipeline] the layered assembly of $OUT stopped before finishing (no $OUT.final.rc): rerun to redo assemble" >&2; exit 1
    fi
    layered_report; rm -f "$OUT.final.pid"   # the pid file goes last: a run killed in the post-steps redoes them
  }
  time_record() { { cat "$OUT.$1.time.tmp" 2>/dev/null | grep -E '^(wall_s|peak_rss_kb)' || true; printf 'exit\t%s\n' "$2"; } > "$OUT.$1.time"; rm -f "$OUT.$1.time.tmp"; }
  all_bg_stop() {
    if [ -n "$scan_pid" ]; then kill_tree "$scan_pid"; fi
    # only a layered run that is still going is stopped, and then assemble must run again (its record goes); a finished one
    # keeps its pid file, so the rerun's all_layers_wait still applies the chain score and ptc steps to PREFIX.final.gtf
    if [ -s "$OUT.final.pid" ] && final_alive; then
      kill_tree "$(cut -d' ' -f1 "$OUT.final.pid")"; rm -f "$OUT.final.pid" "$OUT".final_run.* "$OUT.assemble.env"
    fi
  }
  all_guided() {
    guided_due=0
    all_layers_wait   # register 1934: the guided rule reads $OUT.final.gtf, so the layered assembly has to be whole first
    if [ ! -s "$OUT.guided.gtf" ] || [ "$OUT.guided.gtf" -ot "$OUT.gtf" ] || [ "$OUT.guided.gtf" -ot "$GUIDE" ] || { [ -s "$OUT.final.gtf" ] && [ "$OUT.guided.gtf" -ot "$OUT.final.gtf" ]; }; then
      say "all: guided"; bash "$0" guided "${args[@]}" || { all_bg_stop; say "all: guided FAILED; see $OUT.guided.log"; exit 1; }
    fi
  }
  # rows 1721 / 1722, default since 2026-11-14 (user): PREFIX.earlyends.gtf is a side product of assemble like the guided
  # assembly, outside the plan (nothing reads it, so it never makes a stage rerun) and written last so no other product waits
  # for its BAM pass. RUSTLE_EARLYENDS=0 turns it off. A failure is reported and leaves every other product complete.
  # shortcut: runs after the last stage (A119b +168 s at 72 threads); move it beside families if that wall matters
  # --names: named copies of PREFIX.gtf / PREFIX.final.gtf, outside the plan like the side assemblies (nothing reads them)
  all_names() {
    [ -n "${NAMES:-}" ] && [ -s "$OUT.gtf" ] || return 0   # :- as the dispatch tests run stage_all on its own
    if [ ! -s "$OUT.named.gtf" ] || [ "$OUT.named.gtf" -ot "$OUT.gtf" ] || [ "$OUT.named.gtf" -ot "$NAMES" ] || { [ -s "$OUT.final.gtf" ] && [ ! "$OUT.final.named.gtf" -nt "$OUT.final.gtf" ]; }; then
      say "all: names"
      bash "$0" names "${args[@]}" || say "all: names FAILED (see $OUT.names.log); every other product is complete"
    fi
  }
  all_earlyends() {
    [ "${RUSTLE_EARLYENDS:-1}" != 0 ] && [ -s "$OUT.gtf" ] || return 0
    if [ ! -s "$OUT.earlyends.gtf" ] || [ "$OUT.earlyends.gtf" -ot "$OUT.gtf" ]; then
      say "all: earlyends"
      bash "$0" earlyends "${args[@]}" || say "all: earlyends FAILED (see $OUT.earlyends.log); $OUT.earlyends.gtf is not written, every other product is complete"
    fi
  }
  # rows 1717 / 1720, default since 2026-11-14 (user): PREFIX.capstarts.gtf, a side product like earlyends. Human-tested only;
  # the stage's self-check refuses a library without usable cap G (both gorilla libraries), which is reported and not a failure.
  # A refusal newer than PREFIX.gtf is not retried, so a rerun of `all` does not repeat the BAM pass.
  # rows 1766 / 1806, default since 2026-11-20 (user): the one-read tiers as side products like earlyends, outside the plan
  # (nothing reads them): PREFIX.oneread.gtf (strict) and PREFIX.oneread.ranked.gtf (every one-read chain the default lacks,
  # scored; a ranking, not a filter), plus PREFIX.oneread.promoted.gtf with partner libraries. Costs a second assemble at a
  # floor of one read and two BAM passes. RUSTLE_ONEREAD_TIER=0 turns it off. A failure leaves every other product complete.
  all_oneread() {
    [ "${RUSTLE_ONEREAD_TIER:-1}" != 0 ] && [ -s "$OUT.gtf" ] || return 0
    if [ -s "$OUT.oneread.ranked.gtf" ] && [ ! "$OUT.oneread.ranked.gtf" -ot "$OUT.gtf" ]; then return 0; fi
    say "all: oneread"
    bash "$0" oneread "${args[@]}" || say "all: oneread FAILED (see $OUT.oneread.log); every other product is complete"
  }
  all_capstarts() {
    [ "${RUSTLE_CAPSTARTS:-1}" != 0 ] && [ -s "$OUT.gtf" ] || return 0
    if [ -s "$OUT.capstarts.gtf" ] && [ ! "$OUT.capstarts.gtf" -ot "$OUT.gtf" ]; then return 0; fi
    if [ ! "$OUT.capstarts.log" -ot "$OUT.gtf" ] && grep -q '^capstarts refused' "$OUT.capstarts.log" 2>/dev/null; then return 0; fi
    say "all: capstarts"
    bash "$0" capstarts "${args[@]}" || say "all: capstarts wrote no $OUT.capstarts.gtf (see $OUT.capstarts.log); every other product is complete"
  }
  for item in "${plan[@]}"; do
    st=${item%%:*}; prod=${item#*:}; env_now=$(stage_all_env "$st")
    if [ "$guided_due" = 1 ] && [ "$st" != assemble ]; then all_guided; fi
    if [ "$st" = index ] && [ -n "$idx_pid" ]; then
      ran=1; rm -f "$OUT.$st.env"; say "all: index (waiting for the background build)"
      local rc=0; wait "$idx_pid" || rc=$?; cat "$OUT.index.bg.err" >&2; rm -f "$OUT.index.bg.err"
      time_record "$st" "$rc"
      if [ "$rc" != 0 ]; then
        all_bg_stop   # round 41
        say "all: $st FAILED (exit $rc); see $OUT.$st.time and $OUT.index.log"; exit "$rc"
      fi
      printf '%s\n' "$env_now" > "$OUT.$st.env"; continue
    fi
    if [ "$st" = flag ] && [ -n "$scan_pid" ]; then
      local src=0; wait "$scan_pid" || src=$?; cat "$OUT.flag_scan.bg.err" >&2; rm -f "$OUT.flag_scan.bg.err"
      if [ "$flag_bg" = 1 ]; then
        # round 55: the background job was the whole stage (scan, align, verdict); its wall is the align + verdict part
        ran=1; rm -f "$OUT.$st.env"
        time_record "$st" "$src"
        [ "$src" = 0 ] || { all_bg_stop; say "all: flag FAILED in the background (exit $src); see $OUT.flag_scan.log and $OUT.flag.log"; exit "$src"; }
        say "all: flag ran in the background"; printf '%s\n' "$env_now" > "$OUT.$st.env"; continue
      fi
      [ "$src" = 0 ] || { all_bg_stop; say "all: flag scan FAILED in the background (exit $src); see $OUT.flag_scan.log"; exit "$src"; }
      export RP_FLAG_SCAN_DONE=1
    fi
    if [ "$started" = 0 ]; then
      if [ "$st" = "$FROM_STAGE" ]; then started=1; else say "all: $st skipped (--from $FROM_STAGE)"; continue; fi
    fi
    if [ "$RESUME" = 1 ] && [ "$ran" = 0 ] && [ -s "$prod" ] && [ -f "$OUT.$st.env" ] && [ "$(cat "$OUT.$st.env")" = "$env_now" ] \
       && grep -qx $'exit\t0' "$OUT.$st.time" 2>/dev/null && [ ! "$prod" -nt "$OUT.$st.time" ]; then
      say "all: $st SKIPPED ($prod present, settings unchanged)"
      # --to STAGE must stop here too when the stage itself was up to date (before 2026-11-20 a resumed `--to` run went on)
      if [ -n "$TO_STAGE" ] && [ "$st" = "$TO_STAGE" ]; then [ "$guided_due" = 0 ] || all_guided; all_layers_wait; all_names; all_earlyends; all_capstarts; all_oneread; say "all: stopped after $st (--to $TO_STAGE); product $prod"; return 0; fi
      continue
    fi
    ran=1
    rm -f "$OUT.$st.env" "$OUT.all.env"   # all.env: the single record of builds before 2026-11-07
    say "all: $st"
    local rc=0
    if [ -n "$gtime" ]; then
      "$gtime" -f $'wall_s\t%e\npeak_rss_kb\t%M' -o "$OUT.$st.time.tmp" bash "$0" "$st" "${args[@]}" || rc=$?
    else
      # compute nodes without GNU time (r1543 e2e, 2026-11-06): wall only, no peak RSS
      local st0=$SECONDS; bash "$0" "$st" "${args[@]}" || rc=$?; printf 'wall_s\t%s\n' "$(( SECONDS - st0 ))" > "$OUT.$st.time.tmp"
    fi
    time_record "$st" "$rc"
    if [ "$rc" != 0 ]; then
      # round-41 audit: a background scan must not outlive the run (it writes PREFIX.flag_scan.* in place; a rerun's own scan
      # would race it). A background index build is left to finish: it writes a temporary name and renames.
      all_bg_stop
      say "all: $st FAILED (exit $rc); see $OUT.$st.time and the stage's log"; exit "$rc"
    fi
    printf '%s\n' "$env_now" > "$OUT.$st.env"
    # round-39 audit: assemble just ran, so flag will run: its scan starts now, under families and assign. Not under --low-mem
    # (round 40): that profile bounds memory, and the scan (A119b 16.6 GB) would sit beside families' all-vs-all.
    if [ "$st" = assemble ] && [ "$LOW_MEM" = 0 ] && [[ " ${plan[*]} " == *" flag:"* ]] && [ -z "$FROM_STAGE$TO_STAGE" ]; then
      # round 55: align + verdict read the scan, the loci and the index, nothing of families or assign, so the whole stage
      # runs there (A119b: 42 s no longer after assign) unless the verdict waits for something: the candidates table
      # (--candidates) or an index still being built in the background (a genome's first run)
      if [ "$CANDIDATES" = 0 ] && [ -z "$idx_pid" ]; then
        flag_bg=1; bash "$0" flag_bg "${args[@]}" 2> "$OUT.flag_scan.bg.err" & scan_pid=$!
        say "all: flag (scan, align, verdict) running in the background"
      else
        bash "$0" flag_scan "${args[@]}" 2> "$OUT.flag_scan.bg.err" & scan_pid=$!
        say "all: flag scan running in the background"
      fi
    fi
    if [ -n "$TO_STAGE" ] && [ "$st" = "$TO_STAGE" ]; then [ "$guided_due" = 0 ] || all_guided; all_layers_wait; all_names; all_earlyends; all_capstarts; all_oneread; say "all: stopped after $st (--to $TO_STAGE); product $prod"; return 0; fi
  done
  all_layers_wait
  all_names
  all_earlyends
  all_capstarts
  all_oneread
  say "all: done. stage, wall, peak memory:"
  for item in "${plan[@]}"; do
    st=${item%%:*}
    awk -F'\t' -v st="$st" '$1=="wall_s"{w=$2} $1=="peak_rss_kb"{m=$2} END{printf "[rustle_pipeline]   %-11s %8.0f s  %s\n", st, w, m=="" ? "  peak n/a (no GNU time)" : sprintf("%6.1f GB", m/1048576)}' "$OUT.$st.time" >&2
  done
}
case "$STAGE" in
  merged) stage_merged;;
  assemble) stage_assemble; if [ -n "${GUIDE:-}" ] && [ -z "${RP_ALL:-}" ]; then stage_guided; fi; if [ -n "${NAMES:-}" ] && [ -z "${RP_ALL:-}" ]; then stage_names; fi;; guided) stage_guided;; names) stage_names;; eqclasses) stage_eqclasses;; famwise) stage_famwise;; psi) stage_psi;; multicopy) stage_multicopy;; siblings) stage_siblings;; oneread) stage_oneread || exit $?;; flairlike) stage_flairlike;; altends) stage_altends;; trunc5) stage_trunc5;; chainscore) stage_chainscore;; shiftdrop) stage_shiftdrop;; shortexons) stage_shortexons;; capstarts) stage_capstarts;; earlyends) stage_earlyends;; repeatloci) stage_repeatloci;; exonmap) stage_exonmap;; relatives) stage_relatives;; sdcore) stage_sdcore;; famreps) stage_famreps;; families) stage_families;; candidates) stage_candidates;; catalog) stage_catalog;;
  assign) stage_assign;; o2) stage_o2;; o2bench) stage_o2bench;; members) stage_members;; provenance) stage_provenance;; flag) stage_flag;; flag_scan) stage_flag_scan;; index) stage_index;;
  # `all`'s background flag stage: the scan, then align + verdict, whose wall goes to PREFIX.flag.time
  flag_bg) stage_flag_scan; export RP_FLAG_SCAN_DONE=1; fb0=$SECONDS; stage_flag; printf 'wall_s\t%s\n' "$(( SECONDS - fb0 ))" > "$OUT.flag.time.tmp";;
  # --legacy-catalog: `assign` reads the legacy catalog, so `all` builds it; --candidates: `all` runs the (opt-in) candidates stage
  all) stage_all;;
  *) echo "unknown stage $STAGE" >&2; exit 2;;
esac
say "done"
