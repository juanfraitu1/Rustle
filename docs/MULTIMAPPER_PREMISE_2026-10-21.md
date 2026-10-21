# Multi-mapping long reads: what they cost assembly, what O2 can decide, what O3 recovers (2026-10-21)

**The question (user, 2026-10-21).** "I made my entire thesis on the assumption that multi-mapping long reads were really important
information being lost, that minimap2's designation as primary is no better than a coin toss and that stringtie and flair are
incomplete without this important info. Can we actually check without O2 and O3 ...?"

This page collects the answers of registers 1360-1367. Every row was pre-registered, except register 1366, which is descriptive
over existing runs. Every row was checked by the metric-auditor; the numbers below are the audited ones.

## Short answers

1. **Coin-toss primaries are rare.**
   - An exact AS tie between the primary and another alignment happens for 0.14% (gorilla OR6737) / 0.52% (human A119b) of reads.
   - Within 10 AS: 0.54% / 1.30%.
2. **StringTie and FLAIR already behave as primary-only assemblers.**
   - StringTie gets secondaries but never builds a transcript from them: on the lab's slice, a full BAM and a primary-only BAM
     give identical intron chains.
   - FLAIR drops secondaries.
3. **What multi-mapping costs them is small genome-wide, but concentrated in multi-copy loci.**
   - About 3% of expressed annotated transcripts have mostly multi-placed reads, or no primary alignment at all.
   - There the tools recover about half their usual rate, and none of the no-primary ones.
   - Matched for depth and exon count, the shortfall is ~250-340 transcripts per tool per species.
4. **Using secondaries does not recover real copies.**
   - Seeding with secondaries adds ~0.2 points of sensitivity genome-wide; inside the multi-mapping group it adds +34-43%, but at a
     precision cost.
   - O2's certificate on the reads behind those gains: where the read is truly tied, the copy is undecidable for 76% (human) / 83%
     (gorilla) of the gains.
   - Where the read is not tied, minimap2 had already placed it on the other copy.
   - The mainstay genes with a few certified reads (TBC1D3, NPIPB5, SMN2, GOLGA8F, RBMY1) are certified 97-100% to another copy.
5. **The loci where reads multi-map are real duplicated DNA.**
   - Human: 99.1% are in the CHM13 SD90 track (97.6% on ordinary autosomes; 77.3% without satellite SDs).
   - Gorilla: 98.2% are in SEDEF SDs (a post-hoc check).
   - The base rate for genes without tied reads is 15-24%.
   - This is close to true by construction, since reads tie only where sequence is duplicated, and it does not validate any O1
     family.
6. **O2 in a matched genome** (KB3781 fibroblasts on his own T2T):
   - On real reads it is UNDERPOWERED: minimap2 is wrong on only 1-3 of ~107 tied reads that his genome places on a reference copy.
   - With every known copy expressed from both haplotypes it FAILS. Of minimap2's ~500 errors per replicate, the true copy has
     strictly fewer edits on only 3 / 9. The reads do not carry the provenance minimap2 misses.
7. **O3 recovers a copy the reference lacks:**
   - The reference-free chain finds the maternal-only NPIP locus's units, one candidate each, with nothing flagged outside it.
     IsoCon alone gives that copy only 7 of its 37 reads.
   - Once the candidates are added to the reference, alignment gives it 33 of 39 reads (26 of 32 among reads the candidate was not
     built from).
   - Zero false moves among 2,861 other NPIP reads and 951 reads aligned at the candidates' relatives.
   - The O2 ties at that locus drop from 34 to 10.

## What the thesis can claim, and what it should not

| claim | status |
|---|---|
| Multi-mapping is a genome-wide information loss for StringTie / FLAIR | **not supported** (~3% of transcripts, ~2% of the tools' output) |
| In multi-copy loci the tools lose about half their sensitivity | supported (registers 1366, 1367) |
| Adding secondary alignments recovers those copies | **not supported** (gains mostly undecidable or already placed elsewhere) |
| minimap2's primary is a coin toss | true only for exact ties (0.14-0.52% of reads), where the copies are identical over the read and "both" is the right answer |
| O2 recovers provenance the aligner misses | not in a matched genome (registers 1360-1362); its role is certifying placements and abstaining on true ties (registers 1338, 1342) |
| Multi-mapping loci are real duplications, not artefacts | supported for the loci (register 1367); a separate test is needed for the O1 families |
| A copy the reference lacks explains misplaced tied reads, and O3 can name it and return its reads | supported at one locus (registers 1363-1365; not held out) |

## Caveats

- Two testis samples for the assembly analysis (gorilla OR6737, human A119b).
- Human transcript numbers are scored on NCBI RefSeq, not CAT.
- "Matched" means the intron chain matches the annotation; that cannot tell which near-identical copy built it.
- The O2 and O3 tests use one animal (KB3781), with his own assemblies as the answer key only.

## Sources

| register | document |
|---|---|
| 1360, 1361 | `docs/PREREG_o2_corrector_families_2026-10-20.md` (Amendments 1-4, results, audits) |
| 1362 | `docs/PREREG_o2_corrector_idealsim_2026-10-21.md` |
| 1363 | `docs/PREREG_o3_npip_maternal_copy_2026-10-21.md` |
| 1364, 1365 | `docs/PREREG_o3_npip_candidates_realign_2026-10-21.md` |
| 1366 | descriptive, `docs/NEGATIVE_RESULTS_REGISTER.md` (existing runs in `benchmark_collapse/`, `figures/data/fig3_ref_tie.tsv`) |
| 1367 | `docs/PREREG_multicopy_gains_real_2026-10-21.md` |
