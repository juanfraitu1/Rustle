# PREREG ABSENT-COPY DETECTOR ON A DRAFT GENOME (register 1539, 2026-11-05, before any number)

Question (user 2026-11-05 "lets try out the absent on draft regime"): r1537 found that on a draft genome (GGO reads on Kamilah)
every scorable missing copy separates on identity: its reads form a divergence pile at the absorbing locus or stay unmapped,
while present copies' reads align at T2T identity. Can a BLIND rule (draft BAM only, no T2T, no labels) turn that into calls?
Data: `/anvme/workspace/v138da13-rustle/nont2t/gorilla/{nt.bam,read_origin.tsv,truth2.tsv}` (panel reads on Kamilah) and the
human pair as negative control (`.../human/`, GRCh38: the collapsed regime, where the rule should stay silent).
Rule (blind): primary records of nt.bam clustered by landing (same contig, gaps <= 10 kb merged). A DIVERGENCE PILE is a cluster
with >= 10 reads at de > TAU that are >= 25% of the cluster's primary reads; TAU = 0.005, fixed now (HiFi read error ~0.001 plus
a QV-40 draft ~0.0001; the r1537 present gorilla copies sat at de gap <= 0.001 with one at 0.0075, so one present false alarm is
expected). UNPLACED = unmapped reads; reported as a count, not clustered (no sequence clustering in this row).
Calls = divergence piles. Scoring (labels only here): a missing copy is FOUND when >= 50% of its origin reads are high-de reads
of some pile or unmapped; a pile is TRUE when >= 50% of its high-de reads originate from missing copies, else a FALSE ALARM.
Scorable missing copies = >= 10 origin reads (gorilla 17 of 37).
Bars: recall >= 0.5 of scorable missing copies AND precision (true piles / piles) >= 0.5 on gorilla; human control: piles on
present copies' loci <= 2 (the rule must not fire on collapsed copies at TAU, else TAU is a draft-quality artefact).
Stop rule: PASS -> the rule is the registered absent-copy detector for draft genomes (a `flag` input, not a change to flag);
FAIL on recall -> absent copies are visible only with labels, close; FAIL on precision -> list the false piles' origin.
Script: `bench/nont2t/draft_piles.py`.

## Amendment 1 (2026-11-05, before any number)
1. TAU is on absolute de but was justified from the r1537 de GAP; absolute de adds the read's own HiFi error (~0.001-0.002).
   Predicted ceiling: the four missing copies with gap <= 0.0033 (LOC134756368, NPIPB11, LOC129523574, LRPAP1) are expected to
   sit under TAU and be missed, so about 13 of 17 at best. Human control: present GOLGA8 loci sit at gap 0.002-0.003 on GRCh38
   (another individual's haplotype), absolute de ~0.004-0.005, so the control may fail from haplotype divergence; if it does the
   finding is "TAU is not transferable between genome pairs", not a collapsed-copy signal.
2. Clusters are built blind from every primary record of nt.bam (no labels); labels enter only in scoring.
3. The 392 gorilla ambiguous-origin reads (1,751 human) enter the blind rule and are excluded from the TRUE / FALSE and FOUND
   denominators.
4. nt.bam holds panel reads only: the precision here is an upper bound on a whole-sample BAM.
5. No second draft held-out: the KB3781 maternal-reference alignments are not on disk (only the O2 tables), and this row does
   not realign.
6. A pass means "separable at one fixed constant", not a calibrated detector.
